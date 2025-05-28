package projects.gemorna;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.parameters.AbstractSelectionParameter.InconsistentCollectionException;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;

public class MergeGeMoMaGeMoRNA implements JstacsTool {

	public enum Mode{
		UNION,
		INTERSECT,
		INTERMEDIATE,
		ANNOTATE
	}
	
	public enum Confidence{
		HIGH,
		MEDIUM,
		LOW
	}
	
	public static class GFFEntry{
		
		protected String chr;
		protected String source;
		protected String type;
		protected int start;
		protected int end;
		protected double score;
		protected char strand;
		protected int phase;
		protected HashMap<String,String> attributes;
		
		public GFFEntry(String line) {
			String[] parts = line.split("\t");
			this.chr = parts[0];
			this.source = parts[1];
			this.type = parts[2];
			this.start = Integer.parseInt(parts[3]);
			this.end = Integer.parseInt(parts[4]);
			this.score = ".".equals(parts[5]) ? Double.NaN : Double.parseDouble(parts[5]);
			this.strand = parts[6].charAt(0);
			this.phase = ".".equals(parts[7]) ? -1 : Integer.parseInt(parts[7]);
			String[] attrs = parts[8].split(";");
			this.attributes = new HashMap<>();
			for(int i=0;i<attrs.length;i++) {
				String[] keyval = attrs[i].split("=");
				attributes.put(keyval[0].trim(), keyval[1].trim());
				
			}
		}
		
		public GFFEntry(String chr, String source, String type, int start, int end, double score, char strand,
				int phase, HashMap<String, String> attributes) {
			this.chr = chr;
			this.source = source;
			this.type = type;
			this.start = start;
			this.end = end;
			this.score = score;
			this.strand = strand;
			this.phase = phase;
			this.attributes = attributes;
		}
		
		public String toString() {
			
			String attrs = "";
			if(attributes.containsKey("ID")) {
				attrs += "ID"+"="+attributes.get("ID")+";";
			}
			if(attributes.containsKey("Parent")) {
				attrs += "Parent"+"="+attributes.get("Parent")+";";
			}
			for(String key : attributes.keySet()) {
				if(!"ID".equals(key) && !"Parent".equals(key)) {
					attrs += key+"="+attributes.get(key)+";";
				}
			}
			
			return chr+"\t"+source+"\t"+type+"\t"+start+"\t"+end+"\t"+( Double.isNaN(score) ? "." : score )+"\t"+strand+"\t"+( phase<0 ? "." : phase )+"\t"+attrs;
		}
		
		public String getParent() {
			return attributes.get("Parent");
		}
		
		public boolean equalsStartEndStrand(GFFEntry second) {
			return this.start == second.start && this.end == second.end && this.strand == second.strand;
		}
		
	}
	
	public static class GFFComparator implements Comparator<GFFEntry>{
		
		private static GFFComparator instance = new GFFComparator();
		
		private GFFComparator() {
			
		}
		
		@Override
		public int compare(GFFEntry o1, GFFEntry o2) {
			int cmp = Integer.compare(o1.start, o2.start);
			if(cmp == 0) {
				cmp = Integer.compare(o1.end, o2.end);
			}
			return cmp;
		}

		public static Comparator<GFFEntry> getInstance() {
			return instance;
		}
	}
	
	public static class mRNAComparator implements Comparator<mRNA>{
		
		private static mRNAComparator instance = new mRNAComparator();
		
		private mRNAComparator() {
			
		}
		
		@Override
		public int compare(mRNA o1, mRNA o2) {
			int cmp = Integer.compare(o1.cdsStart, o2.cdsStart);
			if(cmp == 0) {
				cmp = Integer.compare(o1.cdsEnd, o2.cdsEnd);
			}
			if(cmp == 0) {
				cmp = o1.cdss.size() - o2.cdss.size();		
			}
			if(cmp == 0) {
				for(int i=0;i<o1.cdss.size() && cmp == 0;i++) {
					cmp = Integer.compare(o1.cdss.get(i).start, o2.cdss.get(i).start);
					if(cmp == 0) {
						cmp = Integer.compare(o1.cdss.get(i).end, o2.cdss.get(i).end);
					}
					
				}
			}
			return cmp;
		}

		public static mRNAComparator getInstance() {
			return instance;
		}
	}
	
	
	public static class Exon extends GFFEntry{
		
		public Exon(GFFEntry entry) {
			super(entry.chr,entry.source,entry.type,entry.start,entry.end,entry.score,entry.strand,entry.phase,entry.attributes);
			
		}
		
		
	}
	
	public static class CDS extends GFFEntry{
		
		public CDS(GFFEntry entry) {
			super(entry.chr,entry.source,entry.type,entry.start,entry.end,entry.score,entry.strand,entry.phase,entry.attributes);
		}
		
	}
	
	public static class mRNA extends GFFEntry{
		
		private ArrayList<Exon> exons;
		private ArrayList<CDS> cdss;
		private ArrayList<GFFEntry> other;
		private int cdsStart;
		private int cdsEnd;
		private boolean strict;
		
		private Confidence conf;
		
		public mRNA(GFFEntry entry,String genePrefix) {
			super(entry.chr,entry.source,entry.type,entry.start,entry.end,entry.score,entry.strand,entry.phase,entry.attributes);
			this.exons = new ArrayList<>();
			this.cdss = new ArrayList<>();
			this.other = new ArrayList<>();
			this.attributes.put("Parent", genePrefix+this.attributes.get("Parent"));
			this.cdsStart = this.cdsEnd = -1;
			this.strict = false;
			this.conf = null;
		}
		
		public void setConfidence(Confidence confidence) {
			this.conf = confidence;
			if(confidence != null) {
				this.attributes.put("confidence", this.conf.toString().toLowerCase());
			}else {
				this.attributes.remove("confidence");
			}
		}
		
		public String getID() {
			return attributes.get("ID");
		}
		
		public boolean isStrict() {
			return strict;
		}
		
		public void setStrict(boolean strict) {
			this.strict = strict;
		}
		
		public void addExon(Exon exon) {
			this.exons.add(exon);
		}
		
		public void addCDS(CDS cds) {
			this.cdss.add(cds);
			if(cdsStart == -1) {
				cdsStart = cds.start;
				cdsEnd = cds.end;
			}else {
				if(cds.start < cdsStart) {
					cdsStart = cds.start;
				}
				if(cds.end>cdsEnd) {
					cdsEnd = cds.end;
				}
			}
		}
		
		public void addOther(GFFEntry other) {
			this.other.add(other);
		}
		

		public void sortInternally() {
			exons.sort(GFFComparator.getInstance());
			cdss.sort(GFFComparator.getInstance());
		}
		
		public boolean equalsCDSs(mRNA second) {
			if(this.cdss.size() != second.cdss.size()) {
				return false;
			}
			for(int i=0;i<cdss.size();i++) {
				if(! cdss.get(i).equalsStartEndStrand(second.cdss.get(i))){
					return false;
				}
			}
			return true;
		}
		
		public boolean equalsExons(mRNA second) {
			if(this.exons.size() != second.exons.size()) {
				return false;
			}
			for(int i=0;i<exons.size();i++) {
				if(! exons.get(i).equalsStartEndStrand(second.exons.get(i))){
					return false;
				}
			}
			return true;
		}

		public boolean equalsCDSBorders(mRNA rnaRNA) {
			return cdsStart == rnaRNA.cdsStart && cdsEnd == rnaRNA.cdsEnd;
		}

		public boolean overlaps(mRNA rnaRNA) {
			return ( start >=rnaRNA.start && start < rnaRNA.end ) || ( rnaRNA.start >= start && rnaRNA.start < end );
		}

		public void print(PrintWriter wr) {
			wr.println(this);
			for(Exon exon : exons) {
				wr.println(exon);
			}
			for(CDS cds : cdss) {
				wr.println(cds);
			}
			for(GFFEntry en : other) {
				wr.println(en);
			}
			
		}
		
	}
	
	public static class Gene extends GFFEntry{
		
		private LinkedList<mRNA> mrnas;
		private LinkedList<GFFEntry> other;

		
		public Gene(GFFEntry entry,String genePrefix) {
			super(entry.chr,entry.source,entry.type,entry.start,entry.end,entry.score,entry.strand,entry.phase,entry.attributes);
			this.mrnas = new LinkedList<>();
			this.attributes.put("ID", genePrefix+this.attributes.get("ID"));
		}
		
		public String getID() {
			return attributes.get("ID");
		}
		
		public void addmRNA(mRNA mrna) {
			this.mrnas.add(mrna);
		}

		public void addOther(GFFEntry other) {
			this.other.add(other);
		}

		public void sortInternally() {
			mrnas.sort(mRNAComparator.getInstance());
			for(mRNA mrna : mrnas) {
				mrna.sortInternally();
			}
		}

		public void print(PrintWriter wr) {
			wr.println(this);
			for(mRNA mrna : mrnas) {
				mrna.print(wr);
			}
		}
		
	}
	
	
	public static HashMap<String,Gene> parseGFF(String filename,String genePrefix) throws IOException {
		
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		
		HashMap<String,Gene> genes = new HashMap<>();
		HashMap<String,mRNA> mrnas = new HashMap<>();
		ArrayList<Exon> exons = new ArrayList<>();
		ArrayList<CDS> cdss = new ArrayList<>();
		ArrayList<GFFEntry> other = new ArrayList<>();
		
		String str = null;
		while( (str = reader.readLine()) != null ) {
			if(str.trim().length()>0 && !str.startsWith("#")) {
				GFFEntry entry = new GFFEntry(str);
				if("exon".equals(entry.type)) {
					exons.add(new Exon(entry));
				}else if("CDS".equals(entry.type)) {
					cdss.add(new CDS(entry));
				}else if("mRNA".equals(entry.type)) {
					mRNA temp = new mRNA(entry,genePrefix);
					mrnas.put(temp.getID(),temp);
				}else if("gene".equals(entry.type)) {
					Gene temp = new Gene(entry,genePrefix);
					genes.put(temp.getID(),temp);
				}else {
					other.add(entry);
				}
			}
		}
		reader.close();
		
		
		for(Exon exon : exons) {
			String parent = exon.getParent();
			if(mrnas.containsKey(parent)) {
				mrnas.get(parent).addExon(exon);
			}else {
				throw new RuntimeException();
			}
		}
		
		for(CDS cds : cdss) {
			String parent = cds.getParent();
			if(mrnas.containsKey(parent)) {
				mrnas.get(parent).addCDS(cds);
			}else {
				throw new RuntimeException();
			}
		}
		
		for(mRNA mrna : mrnas.values()) {
			String parent = mrna.getParent();
			if(genes.containsKey(parent)) {
				Gene temp = genes.get(parent); 
				temp.addmRNA(mrna);
			}else {
				throw new RuntimeException();
			}
		}
		
		for(GFFEntry entry : other) {
			String parent = entry.getParent();
			if(parent != null) {
				if(mrnas.containsKey(parent)) {
					mrnas.get(parent).addOther(entry);
				}else if(genes.containsKey(parent)) {
					genes.get(parent).addOther(entry);
				}
			}
		}
		
		return genes;
		
	}
	
	private static class SplitAndSortedGFF{
		
		private HashMap<String,ArrayList<mRNA>> map;
		private HashMap<String,Gene> genes;
		
		public SplitAndSortedGFF(HashMap<String,Gene> genes) {
			this.genes = genes;
			map = new HashMap<>();
			for(Gene gene : genes.values()) {
				String chr = gene.chr+"#"+gene.strand;
				if(! map.containsKey(chr)) {
					map.put(chr, new ArrayList<mRNA>());
				}
				gene.sortInternally();
				map.get(chr).addAll(gene.mrnas);
			}
			
			for(ArrayList<mRNA> list : map.values()) {
				list.sort(mRNAComparator.getInstance());
			}
			
		}
		
		public SplitAndSortedGFF(HashMap<String,ArrayList<mRNA>> map, HashMap<String,String> geneMap) {
			this.map = map;
			this.genes = new HashMap<>();
			for(String chr : map.keySet()) {
				ArrayList<mRNA> curr = map.get(chr);
				HashMap<String,LinkedList<mRNA>> newGenes = new HashMap<>();
				for(mRNA mrna : curr) {
					String geneID = mrna.getParent();
					if(geneMap.containsKey(geneID)) {
						geneID = geneMap.get(geneID);
						mrna.attributes.put("Parent", geneID);
					}
					if(!newGenes.containsKey(geneID)) {
						newGenes.put(geneID,new LinkedList<mRNA>());
					}
					newGenes.get(geneID).add(mrna);
				}
				for(String key : newGenes.keySet()) {
					LinkedList<mRNA> rnas = newGenes.get(key);
					Gene temp = createGene(rnas,key);
					genes.put(temp.getID(), temp);
				}
			}
		}

		private Gene createGene(LinkedList<mRNA> rnas, String key) {
			String chr = rnas.getFirst().chr;
			String source = rnas.getFirst().source;
			char strand = rnas.getFirst().strand;
			int start = rnas.getFirst().start;
			int end = rnas.getFirst().end;
			
			double maxTie = -1;
			
			for(mRNA mrna : rnas) {
				if(!chr.equals(mrna.chr) || strand != mrna.strand) {
					throw new RuntimeException();//TODO
				}
				if(!source.equals(mrna.source)) {
					source = "merged";
				}
				if(start > mrna.start) {
					start = mrna.start;
				}
				if(end < mrna.end) {
					end = mrna.end;
				}
				if(mrna.attributes.containsKey("tie")) {//TODO further?
					String tie = mrna.attributes.get("tie");
					if(!"NA".equals(tie)) {
						double temp = Double.parseDouble(tie);
						if(temp > maxTie) {
							maxTie = temp;
						}
					}
				}
			}
			
			Gene gene = new Gene(new GFFEntry(chr,source,"gene",start,end,Double.NaN,strand,-1,new HashMap<>()),"");
			gene.attributes.put("transcripts", rnas.size()+"");
			gene.attributes.put("maxTie", (maxTie>=0 ? maxTie : "NA")+"");
			gene.attributes.put("ID", key);
			
			gene.mrnas.addAll(rnas);
			
			return gene;
		}
		
		
		public void print(PrintWriter wr) {
			
			HashMap<String,ArrayList<Gene>> byChr = new HashMap<>();
			for(Gene gene : genes.values()) {
				if(!byChr.containsKey(gene.chr)) {
					byChr.put(gene.chr, new ArrayList<>());
				}
				byChr.get(gene.chr).add(gene);
			}
			
			for(String chr : byChr.keySet()) {
				ArrayList<Gene> genes = byChr.get(chr);
				genes.sort(GFFComparator.getInstance());
				for(Gene gene : genes) {
					gene.sortInternally();
					gene.print(wr);
				}
			}
			
		}
		
	}
	
	
	public static void main(String[] args) throws Exception {
		
		CLI cli = new CLI(new MergeGeMoMaGeMoRNA());
		cli.run(args);
		
	}

	@Override
	public ToolParameterSet getToolParameters() {
		LinkedList<Parameter> pars = new LinkedList<>();
		
		pars.add(new FileParameter("GeMoMa", "GeMoMa predictions", "gff,gff3", true));
		pars.add(new FileParameter("GeMoRNA","GeMoRNA predictions","gff,gff3",true));
		
		try {
			pars.add(new SelectionParameter(DataType.PARAMETERSET, new String[] {"intersect","union","intermediate","annotate"}, new ParameterSet[] {
					new SimpleParameterSet(),
					new SimpleParameterSet(),
					new SimpleParameterSet(
							new FileParameter("GeMoMa-strict","GeMoMa predictions with strict settings","gff,gff3",true),
							new FileParameter("GeMoRNA-strict","GeMoRNA predictions with strict settings","gff,gff3",true)
							),
					new SimpleParameterSet(
							new FileParameter("GeMoMa-strict","GeMoMa predictions with strict settings","gff,gff3",true),
							new FileParameter("GeMoRNA-strict","GeMoRNA predictions with strict settings","gff,gff3",true),
							new SimpleParameter(DataType.BOOLEAN, "Low-confidence", "include low-confidence predictions", true, true)
			)}, new String[] {"Intersection of GeMoRNA and GeMoMa predictions","Union of GeMoRNA and GeMoMa predictions","Intersection and strict GeMoRNA and GeMoMa predictions","Union of GeMoRNA and GeMoMa predictions with annotation as high, medium or low confidence"}, "Mode", "", true));
		} catch (InconsistentCollectionException | IllegalValueException | DatatypeNotValidException e) {
			e.printStackTrace();
		}
		
		
		return new ToolParameterSet(getToolName(), pars);
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		String gemomaFile = (String) parameters.getParameterAt(0).getValue();
		String gemornaFile = (String) parameters.getParameterAt(1).getValue();
		
		SplitAndSortedGFF gemoma = new SplitAndSortedGFF(parseGFF(gemomaFile,"GeMoMa_"));
		SplitAndSortedGFF gemorna = new SplitAndSortedGFF(parseGFF(gemornaFile,"GeMoRNA_"));
		
		int modei = ((SelectionParameter)parameters.getParameterAt(2)).getSelected();
		
		SplitAndSortedGFF gemomaStrict = null;
		SplitAndSortedGFF gemornaStrict = null;
		
		boolean includeLow = true;
		
		if(modei == 2||modei==3) {
			ParameterSet ps= (ParameterSet) parameters.getParameterAt(2).getValue();
			gemomaStrict = new SplitAndSortedGFF(parseGFF( (String)ps.getParameterAt(0).getValue() ,"GeMoMa_"));
			gemornaStrict = new SplitAndSortedGFF(parseGFF( (String)ps.getParameterAt(1).getValue() ,"GeMoRNA_"));
			
			annotateStrict(gemoma,gemomaStrict);
			annotateStrict(gemorna,gemornaStrict);
			
			if(modei == 3) {
				includeLow = (boolean) ps.getParameterAt(2).getValue();
			}
			
		}
		
		
		
		
		
		LinkedList<Result> ress = new LinkedList<>();
		
		Mode mode = null;
		switch(modei) {
			case 0: mode = Mode.INTERSECT; break;
			case 1: mode = Mode.UNION; break;
			case 2: mode = Mode.INTERMEDIATE; break;
			default: mode = Mode.ANNOTATE;
		}
		
		
		SplitAndSortedGFF union = merge(gemoma,gemorna, mode,includeLow);
		
		File temp = File.createTempFile("merge", "gff", new File("."));
		temp.deleteOnExit();
		PrintWriter wr = new PrintWriter(temp);
		union.print(wr);
		wr.close();
		
		TextResult tr = new TextResult("Merged Predictions", "", new FileRepresentation(temp.getAbsolutePath()), true, "gff", getToolName(), "gff", true);
		ress.add(tr);
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(ress), parameters, getToolName(), new Date(System.currentTimeMillis()) );
		
	}

	private void annotateStrict(SplitAndSortedGFF gemoma, SplitAndSortedGFF gemomaStrict) {
		mRNAComparator cmp = mRNAComparator.getInstance();
		
		for(String chr : gemoma.map.keySet()) {
			ArrayList<mRNA> ma = gemoma.map.get(chr);
			ArrayList<mRNA> strict = gemomaStrict.map.containsKey(chr) ? gemomaStrict.map.get(chr) : new ArrayList<>();
			
			int i=0;
			int j=0;
			while(i<ma.size() && j<strict.size()) {
				mRNA maRNA = ma.get(i);
				mRNA strictRNA = strict.get(j);
				
				int cmpRes = cmp.compare(maRNA, strictRNA);
				
				if( cmpRes == 0) {
					maRNA.setStrict(true);
					i++;
					j++;
				}else { 
					if(cmpRes < 0) {
						i++;
					}else if(cmpRes > 0) {
						j++;
					}
				}
			}
		}
		
		
	}

	private SplitAndSortedGFF merge(SplitAndSortedGFF gemoma, SplitAndSortedGFF gemorna, Mode mode, boolean includeLow) {
		
		HashMap<String,String> geneMap = new HashMap<>();
		
		HashMap<String,ArrayList<mRNA>> joinedMap = new HashMap<>();
		
		mRNAComparator cmp = mRNAComparator.getInstance();
		
		HashSet<String> keys = new HashSet(gemoma.map.keySet());
		keys.addAll(gemorna.map.keySet());
		
		for(String chr : keys) {
			ArrayList<mRNA> ma = gemoma.map.containsKey(chr) ? gemoma.map.get(chr) : new ArrayList<>();
			ArrayList<mRNA> rna = gemorna.map.containsKey(chr) ? gemorna.map.get(chr) : new ArrayList<>();
			ArrayList<mRNA> joinedList = new ArrayList<>();
			int i=0;
			int j=0;
			while(i<ma.size() && j<rna.size()) {
				mRNA maRNA = ma.get(i);
				mRNA rnaRNA = rna.get(j);
				
				int cmpRes = cmp.compare(maRNA, rnaRNA);
				
				if( cmpRes == 0) {
					mRNA joined = join(maRNA,rnaRNA,geneMap);
					if(mode == Mode.ANNOTATE) {
						joined.setConfidence(Confidence.HIGH);
					}
					joinedList.add(joined);
					i++;
					j++;
				}else { 
					if(maRNA.overlaps(rnaRNA)) {
						addToMap(maRNA,rnaRNA,geneMap);
					}
					if(cmpRes < 0) {
						if(mode == Mode.UNION || (mode == Mode.ANNOTATE && includeLow) || maRNA.isStrict()) {
							if(mode == Mode.ANNOTATE) {
								maRNA.setConfidence(maRNA.isStrict() ? Confidence.MEDIUM : Confidence.LOW);
							}
							joinedList.add(maRNA);
						}
						i++;
					}else if(cmpRes > 0) {
						if(mode == Mode.UNION || (mode == Mode.ANNOTATE && includeLow) || rnaRNA.isStrict()) {
							if(mode == Mode.ANNOTATE) {
								rnaRNA.setConfidence(rnaRNA.isStrict() ? Confidence.MEDIUM : Confidence.LOW);
							}
							joinedList.add(rnaRNA);
						}
						j++;
					}
				}
			}
			while(i<ma.size()) {
				mRNA maRNA = ma.get(i);
				if(mode == Mode.UNION || (mode == Mode.ANNOTATE && includeLow) || maRNA.isStrict()) {
					if(mode == Mode.ANNOTATE) {
						maRNA.setConfidence(maRNA.isStrict() ? Confidence.MEDIUM : Confidence.LOW);
					}
					joinedList.add(maRNA);
				}
				i++;
			}
			while(j<rna.size()) {
				mRNA rnaRNA = rna.get(j);
				if(mode == Mode.UNION || (mode == Mode.ANNOTATE && includeLow) || rnaRNA.isStrict()) {
					if(mode == Mode.ANNOTATE) {
						rnaRNA.setConfidence(rnaRNA.isStrict() ? Confidence.MEDIUM : Confidence.LOW);
					}
					joinedList.add(rnaRNA);
				}
				j++;
			}
			joinedMap.put(chr, joinedList);
		}
		
		return new SplitAndSortedGFF(joinedMap,geneMap);
		
	}

	private void addToMap(mRNA maRNA, mRNA rnaRNA, HashMap<String, String> geneMap) {
		geneMap.put(rnaRNA.getParent(), maRNA.getParent());		
	}

	private mRNA join(mRNA maRNA, mRNA rnaRNA, HashMap<String, String> geneMap) {
		mRNA joined = new mRNA(rnaRNA,"");
		joined.source = "merged";
		joined.attributes = new HashMap<>();
		
		joined.attributes.put("ID", maRNA.getID());
		joined.attributes.put("Parent", maRNA.getParent());
		
		String altStr = rnaRNA.getID();
		if(rnaRNA.attributes.containsKey("alternative")) {
			String temp = rnaRNA.attributes.get("alternative");
			temp = temp.replaceAll("\"", "");
			altStr += ","+temp;
		}
		if(maRNA.attributes.containsKey("alternative")) {
			String temp = maRNA.attributes.get("alternative");
			temp = temp.replaceAll("\"", "");
			altStr += ","+temp;
		}
		joined.attributes.put("alternative", "\""+altStr+"\"");
		
		for(String key : maRNA.attributes.keySet()) {
			if(! ("alternative".equals(key) || "ID".equals(key) || "Parent".equals(key) ) ){
				joined.attributes.put("GeMoMa_"+key, maRNA.attributes.get(key));
			}
		}
		for(String key : rnaRNA.attributes.keySet()) {
			if(! ("alternative".equals(key) || "ID".equals(key) || "Parent".equals(key) ) ){
				joined.attributes.put("GeMoRNA_"+key, rnaRNA.attributes.get(key));
			}
		}
		
		geneMap.put(rnaRNA.getParent(), maRNA.getParent());
		
		for(Exon exon : rnaRNA.exons) {
			Exon joinedExon = new Exon(exon);
			joinedExon.attributes.put("Parent", maRNA.getID());
			joined.addExon(joinedExon);
		}
		for(CDS cds : maRNA.cdss) {
			CDS joinedCDS = new CDS(cds);
			joinedCDS.attributes.put("Parent", maRNA.getID());
			joined.addCDS(joinedCDS);
		}
		for(GFFEntry en : rnaRNA.other) {
			GFFEntry joinedOther = new GFFEntry(en.chr, en.source, en.type, en.start, en.end, en.score, en.strand, en.phase, en.attributes);
			joinedOther.attributes.put("Parent", maRNA.getID());
			joined.addOther(joinedOther);
		}
		
		return joined;
		
	}

	@Override
	public String getToolName() {
		return "Merge";
	}

	@Override
	public String getToolVersion() {
		return "1.2";
	}

	@Override
	public String getShortName() {
		return "merge";
	}

	@Override
	public String getDescription() {
		return "";
	}

	@Override
	public String getHelpText() {
		return "";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}

	@Override
	public ToolResult[] getTestCases(String path) {
		return null;
	}

	@Override
	public void clear() {
		
	}

	@Override
	public String[] getReferences() {
		return null;
	}

}
