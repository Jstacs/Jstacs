package projects.gemorna;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Locale;

import de.jstacs.algorithms.optimization.ConstantStartDistance;
import de.jstacs.algorithms.optimization.DifferentiableFunction;
import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.TerminationException;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.io.ArrayHandler;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SafeOutputStream;
import de.jstacs.utils.ToolBox;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import projects.gemoma.ExtractRNAseqEvidence.Stranded;
import projects.gemoma.Tools;
import projects.gemorna.ReadGraph.Edge;
import projects.gemorna.ReadGraph.Node;

public class SplicingGraph {

	
	private static class StopIterator implements Iterator<int[]>{

		private char[] sequence;
		private int offset;
		private int len;
		private int start;
		
		public void set(char[] sequence, int offset) {
			this.sequence = sequence;
			this.offset = offset;
			len = -1;
			start = -3;
		}
		
		public boolean isAtEnd() {
			return len < 0 && offset+3 > sequence.length;
		}
		
		@Override
		public boolean hasNext() {
			if(len >= 0) {
				return true;
			}
			int localOff = offset;
			if(offset+3 > sequence.length) {
				return false;
			}
			for( ; offset+3 <= sequence.length; offset+=3 ) {
				if(Genome.isStop[sequence[offset]][sequence[offset+1]][sequence[offset+2]]) {
					len = offset-localOff;
					offset += 3;
					return true;
				}else if(Genome.isStart[sequence[offset]][sequence[offset+1]][sequence[offset+2]]) {
					if(start < 0) {
						start = offset-localOff;
					}
				}
				
			}
			len = offset-localOff;
			return true;
		}

		@Override
		public int[] next() {
			int[] temp = new int[] {len/3,start/3};
			len = -1;
			start = -3;
			return temp;
		}
		
		
		
		
	}
	
	
	private static class CDSCandidateList{
		
		private int[] lens;
		private int[] frames;
		private int[] idxs;
		private double[] vals;
		private boolean[] lasts;
		
		
		private CDSCandidateList(int numBest) {
			lens = new int[numBest];
			frames = new int[numBest];
			idxs = new int[numBest];
			vals = new double[numBest];
			lasts = new boolean[numBest];
		}
		
		
		private void add(int len, int frame, int numCDSIntrons, int numIntrons, int idx, boolean last) {
			double val = getValue(len,frame,numCDSIntrons,numIntrons,idx);
			int i = vals.length;
			while(i>0 && val > vals[i-1]) {
				i--;
			}
			if(i<vals.length) {
				for(int j=i;j<vals.length-1;j++) {
					vals[j+1] = vals[j];
					lens[j+1] = lens[j];
					frames[j+1] = frames[j];
					idxs[j+1] = idxs[j];
					lasts[j+1] = lasts[j];
				}
				vals[i] = val;
				lens[i] = len;
				frames[i] = frame;
				idxs[i] = idx;
				lasts[i] = last;
			}
		}
		
		
		private double getValue(int len, int frame, int numCDSIntrons, int numIntrons, int idx) {
			return len*(1.0 + 0.1*numCDSIntrons + 0.1*(numCDSIntrons+1.0)/(numIntrons+1.0));
		}


		public int[] getBest() {
			if(vals[0]>0) {
				return new int[] {lens[0],frames[0],idxs[0],lasts[0] ? 1 : 0};
			}else {
				return new int[] {0,0,0,-1};
			}
		}
		
		
	}
	
	private static class IntronKey{
		
		private static int hashInt = (int)'#';
		
		private int start;
		private int end;
		
		public IntronKey(int start, int end) {
			this.start = start;
			this.end = end;
		}
		
		public void setStartEnd(int start, int end) {
			this.start = start;
			this.end = end;
		}
		
		@Override
		public boolean equals(Object o) {
			return ((IntronKey)o).start == start && ((IntronKey)o).end == end;
		}
		
		@Override
		public int hashCode() {
			int value = start;
			value = value * 31 + hashInt;
			value = value * 31 + end;
			return value;
			
		}
		
	}
	
	private static NumberFormat df = DecimalFormat.getInstance(Locale.ENGLISH);
	static {
		df.setMaximumFractionDigits(3);
		df.setGroupingUsed(false);
	}
	
	private static int SHORTEXON = 10;
	
	ContiguousRegion[] contRegs;
	LinkedList<ComparableElement<Integer,Integer>>[] adjList;
	private double[][] adjCounts;
	
	Intron[] introns;
	
	private String chromosome;
	private int regionStart;
	
	private HashMap<Integer,Character> strands;
	
	
	private static int[] intersectSortedArrays(int[] edge, int[] left) {
		int i=0,j=0;
		IntList res = new IntList();
		while(i < edge.length && j < left.length) {
			if(left[j] == edge[i]) {
				res.add(left[j]);
				i++;
				j++;
			}else if(left[j] < edge[i]) {
				j++;
			}else if(left[j] > edge[i]) {
				i++;
			}
		}
		return res.toArray();
	}
	
	private static void unionSortedArrays(int[] res, int[] edge, int[] left) {
		int i=0,j=0;
		int k=0;
		while(i < edge.length && edge[i] > -1 && j < left.length && left[j] > -1) {
			if(left[j] == edge[i]) {
				res[k] = left[j];
				k++;
				i++;
				j++;
			}else if(left[j] < edge[i]) {
				res[k]=left[j];
				k++;
				j++;
			}else if(left[j] > edge[i]) {
				res[k] = edge[i];
				k++;
				i++;
			}
		}
		while(i<edge.length && edge[i] > -1) {
			res[k] = edge[i];
			k++;
			i++;
		}
		while(j<left.length && left[j] > -1) {
			res[k] = left[j];
			k++;
			j++;
		}
		res[k] = -1;
	}
	
	private static int[] setdiff(int[] ref, int[] diff) {
		int i=0,j=0;
		IntList res = new IntList();
		while(i < ref.length && ref[i] > -1 && j < diff.length && diff[j] > -1) {
			if(diff[j] == ref[i]) {
				i++;
				j++;
			}else if(diff[j] < ref[i]) {
				j++;
			}else if(diff[j] > ref[i]) {
				res.add(ref[i]);
				i++;
			}
		}
		while(i<ref.length && ref[i] > -1) {
			res.add(ref[i]);
			i++;
		}
		return res.toArray();
	}
	
	
	public SplicingGraph(ReadGraph graph, Region reg, int minIntronLength, int[] proposedSplits) throws CloneNotSupportedException {
		this.chromosome = graph.chrom;
		this.regionStart = graph.regionStart;
		
		
		LinkedList<SplicingGraph.ContiguousRegion> contRegs = new LinkedList<SplicingGraph.ContiguousRegion>();
		LinkedList<SplicingGraph.Intron> introns = new LinkedList<SplicingGraph.Intron>();
		
		ContiguousRegion[] regs = new ContiguousRegion[graph.nodes.length];
		
		IntList currNodes = new IntList();
		
		for(int i=0;i<graph.nodes.length;i++) {
			if(graph.nodes[i].getNumberOfIncomingEdges() == 0 && graph.nodes[i].getNumberOfOutgoingEdges()==0) {
				currNodes.clear();
				continue;
			}
			
			if(currNodes.contains(i)==-1) {
				currNodes.add(i);
			}

			if(graph.nodes[i].getNumberOfOutgoingEdges() == 1 &&
					graph.nodes[i].getOutgoingEdges().iterator().next().getRelEnd() == i+1 &&
					graph.nodes[i+1].getNumberOfIncomingEdges() == 1 
					//!Arrays.stream(proposedSplits).anyMatch(n->n==(j+regionStart)) 
							&& Arrays.binarySearch(proposedSplits, i+1+regionStart)<0
					) {
						
						//System.out.println("true");
						currNodes.add(i+1);
			}else {
				if(currNodes.length()>0) {
					ContiguousRegion region = new ContiguousRegion(currNodes.get(0),currNodes.get(currNodes.length()-1),getNodes(graph.nodes,currNodes));
					contRegs.add(region);
					if(regs[region.regionStart] != null || regs[region.regionEnd] != null) {
						throw new RuntimeException("Overlapping contregs");
					}
					regs[region.regionStart] = regs[region.regionEnd] = region;
				}
				currNodes.clear();
				if(i+1<graph.nodes.length) {
					currNodes.add(i+1);
				}
			}
		}
		if(currNodes.length()>0) {
			ContiguousRegion region = new ContiguousRegion(currNodes.get(0),currNodes.get(currNodes.length()-1),getNodes(graph.nodes,currNodes));
			contRegs.add(region);
		}
		
		
		Iterator<Edge> edgeIt = graph.edges.iterator();
		
		while(edgeIt.hasNext()) {
			Edge e = edgeIt.next();
			int start = e.getRelStart();
			int end = e.getRelEnd();
			if(regs[start] != null && regs[end] != null) {
				if(regs[start] == regs[end]) {
					continue;
				}
				Intron intron = new Intron(regs[start], regs[end], e);
				regs[start].right.add(intron);
				regs[end].left.add(intron);
				introns.add(intron);
			}
		}
		
		LinkedList<ContiguousRegion> crRemove = new LinkedList<SplicingGraph.ContiguousRegion>();
		LinkedList<Intron> inRemove = new LinkedList<SplicingGraph.Intron>();
		
		for(ContiguousRegion region : contRegs) {
			if(region.regionEnd-region.regionStart+1 <= SHORTEXON ) {
				if(region.left.size() == 0 && region.regionEnd+1 < regs.length && regs[region.regionEnd+1] != null && regs[region.regionEnd+1].left.size()>1) {
					crRemove.add(region);
					inRemove.addAll(region.right);
				}
				if(region.right.size() == 0 && region.regionStart-1>=0 && regs[region.regionStart-1] != null && regs[region.regionStart-1].right.size()>1) {
					crRemove.add(region);
					inRemove.addAll(region.left);
				}
			}
			if(region.left.size() == 0 && region.right.size() == 0) {
				int start = region.regionStart;
				int end = region.regionEnd;
				innerloop:
				for(Intron intron : introns) {
					if(intron.left.regionEnd < start && intron.right.regionStart > end) {
						if(intron.edge.getNumberOfReads()*0.5 > region.getNumberOfReads()) {
							crRemove.add(region);
							break innerloop;
						}
					}
				}
			}
		}
		
		for(ContiguousRegion region : crRemove) {
			regs[region.regionStart] = regs[region.regionEnd] = null;
		}
		
		for(ContiguousRegion region: contRegs) {
			region.left.removeAll(inRemove);
			region.right.removeAll(inRemove);
		}
		
		contRegs.removeAll(crRemove);
		introns.removeAll(inRemove);
		
		
		this.contRegs = contRegs.toArray(new ContiguousRegion[0]);
		this.introns = introns.toArray(new Intron[0]);
		

		addReads(reg, minIntronLength);
		
		crRemove.clear();
		inRemove.clear();
		
		for(int i=0;i<this.contRegs.length;i++) {
			if(this.contRegs[i].readIdxs.length == 0) {
				crRemove.add(this.contRegs[i]);
				inRemove.addAll(this.contRegs[i].left);
				inRemove.addAll(this.contRegs[i].right);
			}
		}
		for(int i=0;i<this.introns.length;i++) {
			if(this.introns[i].readIdxs.length == 0 //TODO start modify
					|| (this.introns[i].getLength()>0 && !this.introns[i].isCanonical(minIntronLength))
					//TODO end modify
					) {
				inRemove.add(this.introns[i]);
				for(int j=0;j<this.contRegs.length;j++) {
					this.contRegs[j].left.remove(this.introns[i]);
					this.contRegs[j].right.remove(this.introns[i]);
				}
			}
		}
		
		
		for(ContiguousRegion region : crRemove) {
			regs[region.regionStart] = regs[region.regionEnd] = null;
		}
		
		contRegs.removeAll(crRemove);
		introns.removeAll(inRemove);
		
		this.contRegs = contRegs.toArray(new ContiguousRegion[0]);
		this.introns = introns.toArray(new Intron[0]);
			
			
		this.adjCounts = new double[this.contRegs.length][this.contRegs.length];
		this.adjList = ArrayHandler.createArrayOf(new LinkedList<ComparableElement<Integer,Integer>>(), this.contRegs.length);
		
		for(int i=0;i<this.introns.length;i++) {
			int leftIdx = contRegs.indexOf( this.introns[i].left );
			int rightIdx = contRegs.indexOf(this.introns[i].right);
			this.adjList[leftIdx].add( new ComparableElement<Integer, Integer>(rightIdx, i) );
		}
		for(int i=0;i<this.adjList.length;i++) {
			Collections.sort(this.adjList[i],new Comparator<ComparableElement<Integer,Integer>>() {

				@Override
				public int compare(ComparableElement<Integer,Integer> o1, ComparableElement<Integer,Integer> o2) {
					return Integer.compare(SplicingGraph.this.introns[o2.getWeight()].readIdxs.length,SplicingGraph.this.introns[o1.getWeight()].readIdxs.length);
				}
				
			});
		}
		
		
		
	}	
	
	
	SplicingGraph(projects.gemoma.Analyzer.Transcript t) {
		this.chromosome = t.getChromosome();
		this.regionStart = t.getMin();
		this.contRegs = new ContiguousRegion[t.getNumberOfParts()];
		for(int idx=0;idx<t.getNumberOfParts();idx++) {
			int[] part = t.getPart(idx);
			this.contRegs[idx] = new ContiguousRegion(part[0]-this.regionStart,part[1]-this.regionStart,new LinkedList<Node>());
		}
		Arrays.sort(this.contRegs,new Comparator<ContiguousRegion>() {

			@Override
			public int compare(ContiguousRegion o1, ContiguousRegion o2) {
				return Integer.compare(o1.regionStart, o2.regionStart);
			}
			
		});
		this.introns = new Intron[contRegs.length-1];
		for(int i=1;i<this.contRegs.length;i++) {
			this.introns[i-1] = new Intron(this.contRegs[i-1], this.contRegs[i], new Edge(this.contRegs[i-1].regionEnd,this.contRegs[i].regionStart));
		}
	}

	Transcript createTranscript(projects.gemoma.Analyzer.Transcript t) {
		LinkedList<Integer> exons = new LinkedList<>();
		LinkedList<Integer> introns = new LinkedList<>();
		for(int i=0;i<t.getNumberOfParts();i++) {
			exons.add(i);
			if(i<this.introns.length) {
				introns.add(i);
			}
		}
		Transcript t2 = new Transcript(exons,introns);
		t2.id = t.getID();
		t2.geneID = t.getParent();
		t2.strand = t.getStrand();
		t2.attributes = t.getInfos();
		return t2;
	}
	
	private void addReads(Region reg, int minIntronLength) throws CloneNotSupportedException {
		
		
		
		int len = reg.getRegionEnd()-regionStart+1;
		int[] regions = new int[len];
		IntList[] regReads = ArrayHandler.createArrayOf(new IntList(), contRegs.length);
		IntList[] firstReads = ArrayHandler.createArrayOf(new IntList(), contRegs.length);
		IntList[] lastReads = ArrayHandler.createArrayOf(new IntList(), contRegs.length);
		IntList[] intronReads = ArrayHandler.createArrayOf(new IntList(), introns.length);
		
		HashMap<IntronKey,Integer> intronMap = new HashMap<IntronKey,Integer>();
		for(int i=0;i<introns.length;i++) {
			//String key = introns[i].edge.getRelStart()+"#"+introns[i].edge.getRelEnd();
			IntronKey key = new IntronKey(introns[i].edge.getRelStart(),introns[i].edge.getRelEnd());
			intronMap.put(key, i);
		}
		
		
		Arrays.fill(regions, -1);
		for(int i=0;i<contRegs.length;i++) {
			for(int j=contRegs[i].regionStart;j<=contRegs[i].regionEnd;j++) {
				if(regions[j]>-1) {
					throw new RuntimeException();
				}
				regions[j] = i;
			}
		}
		
		
		LinkedList<SAMRecord> list = reg.getReads();
		HashMap<String,Integer> idMap = reg.getIDMap();
		this.strands = reg.getStrands();
		
		
		Iterator<SAMRecord> it = list.iterator();
		
		while(it.hasNext()) {
			SAMRecord sr = it.next();
			
			if(sr.getAlignmentEnd() < regionStart || sr.getAlignmentStart()>regionStart+contRegs[contRegs.length-1].regionEnd) {
				continue;
			}
			
			String name = sr.getReadName();
			int idx= idMap.get(name);
			
			addRead(sr, idx, minIntronLength, intronReads, intronMap, 
					regions, regReads, firstReads, lastReads);
		}
		
		
		for(int i=0;i<contRegs.length;i++) {
			regReads[i].sortAndMakeUnique();
			firstReads[i].sortAndMakeUnique();
			lastReads[i].sortAndMakeUnique();
			contRegs[i].setSortedReads(regReads[i].toArray(), firstReads[i].toArray(), lastReads[i].toArray());
		}
		
		for(int i=0;i<introns.length;i++) {
			intronReads[i].sortAndMakeUnique();
			introns[i].readIdxs = intronReads[i].toArray();
		}
		
		
	}
	
	
	private void addRead(SAMRecord read, int readIdx, int minIntronLength, IntList[] intronReads, HashMap<IntronKey,Integer> intronMap, 
			int[] regions, IntList[] regReads, IntList[] firstReads, IntList[] lastReads) {
		Iterator<AlignmentBlock> blockIt = read.getAlignmentBlocks().iterator();
		
		int lastRelEnd = -1;
		
		int lastRegion = -1;
		
		IntronKey key = new IntronKey(0,0);
		
		while(blockIt.hasNext()) {
			
			AlignmentBlock block = blockIt.next();

			int relStart = block.getReferenceStart()-regionStart;
			
			if(lastRelEnd > -1 && relStart >= 0) {
				if(relStart-lastRelEnd < minIntronLength) {
					for(int i=lastRelEnd+1;i<=relStart;i++) {
						key.setStartEnd(i-1, i);
						if(intronMap.containsKey(key)) {
							int idx = intronMap.get(key);
							intronReads[idx].add(readIdx);
						}
					}
				}else {
					key.setStartEnd(lastRelEnd, relStart);
					if(intronMap.containsKey(key)) {
						int idx = intronMap.get(key);
						intronReads[idx].add(readIdx);
					}else {
					}
				}
			}
			
			
			
			for(int i=Math.max(0,-relStart);i<block.getLength() && relStart+i < regions.length;i++) {
				
				int regionIdx = regions[relStart+i];
				if(regionIdx > -1) {
					if(regionIdx != lastRegion) {
						
						regReads[regionIdx].add(readIdx);
					}
					if(relStart+i == contRegs[regionIdx].regionStart) {
						firstReads[regionIdx].add(readIdx);
					}
					if(relStart+i == contRegs[regionIdx].regionEnd) {
						lastReads[regionIdx].add(readIdx);
					}
					lastRegion = regionIdx;
				}
				if(i>0) {
					key.setStartEnd(relStart+i-1, relStart+i);
					if(intronMap.containsKey(key)) {
						int idx = intronMap.get(key);
						intronReads[idx].add(readIdx);
					}
				}
			}
			lastRelEnd = relStart+block.getLength()-1;
		}
	}

	public String toString() {
		StringBuffer sb = new StringBuffer();
		for(ContiguousRegion region : contRegs) {
			sb.append(region.toString(regionStart)+"\n");
		}
		for(Intron intron : introns) {
			sb.append(intron.toString(regionStart)+"\n");
		}
		return sb.toString();
	}
	
	public void getIntersections() {
		for(int i=0;i<contRegs.length;i++) {
			for(int j=0;j<contRegs.length;j++) {
				int n = intersectSortedArrays(contRegs[i].readIdxs, contRegs[j].readIdxs).length;
				System.out.print(n+"\t");
			}
			for(int j=0;j<introns.length;j++) {
				int n = intersectSortedArrays(contRegs[i].readIdxs,introns[j].readIdxs).length;
				System.out.print(n+"\t");
			}
			System.out.println();
		}
		for(int i=0;i<introns.length;i++) {
			for(int j=0;j<contRegs.length;j++) {
				int n = intersectSortedArrays(introns[i].readIdxs, contRegs[j].readIdxs).length;
				System.out.print(n+"\t");
			}
			for(int j=0;j<introns.length;j++) {
				int n = intersectSortedArrays(introns[i].readIdxs,introns[j].readIdxs).length;
				System.out.print(n+"\t");
			}
			System.out.println();
		}
	}
	
	private static LinkedList<Node> getNodes(ReadGraph.Node[] nodes, IntList nodeIdxs) {
		LinkedList<Node> li = new LinkedList<ReadGraph.Node>();
		for(int i=0;i<nodeIdxs.length();i++) {
			li.add(nodes[nodeIdxs.get(i)]);
		}
		return li;
	}
	
	public LinkedList<Gene> finalize(LinkedList<Transcript> list, String geneBase, boolean useChrPrefix, int firstID, double minReadsPerGene, int minProteinLength) {
		
		if(list.size() == 0) {
			return new LinkedList<Gene>();
		}
		
		Collections.sort(list);
		
		for(Transcript t : list) {
			t.addStrandAndCDS(true, minProteinLength,false);
		}
		
		LinkedList<Transcript> gene = new LinkedList<SplicingGraph.Transcript>();
		LinkedList<Gene> allGenes = new LinkedList<SplicingGraph.Gene>();
		
		Transcript curr = null;
		
		for(Transcript t : list) {
			if(curr == null) {
				curr = t;
				gene.add(t);
			}else {
				if(t.getStart() <= curr.getEnd()) {
					if(t.getEnd() > curr.getEnd()) {
						curr = t;
					}
					gene.add(t);
				}else {
					
					firstID = addGenes(allGenes,geneBase,useChrPrefix,firstID,gene,minReadsPerGene);
					gene.clear();
					gene.add(t);
					curr = t;
				}
			}
		}
		
		firstID = addGenes(allGenes,geneBase,useChrPrefix,firstID,gene,minReadsPerGene);
		
		return allGenes;
		
	}
	
	private Gene getGene(String geneBase, boolean useChrPrefix,  int firstID, LinkedList<Transcript> gene, char strand, double minReadsPerGene) {
		
		if(useChrPrefix) {
			geneBase += chromosome+".";
		}
		
		LinkedList<Transcript> sub = new LinkedList<SplicingGraph.Transcript>();
		int j=1;
		double reads = 0.0;
		for(Transcript t : gene) {
			if(t.strand == strand) {
				t.setGene(geneBase+firstID);
				t.setId(geneBase+firstID+"."+j+"");
				j++;
				sub.add(t);
				reads += t.getAbundance();
			}
		}
		if(sub.size() == 0 || reads < minReadsPerGene) {
			return null;
		}else {
			return new Gene(geneBase+firstID,sub,strand);
		}
		
	}
	
	private int addGenes(LinkedList<Gene> allGenes, String geneBase, boolean useChrPrefix,  int firstID, LinkedList<Transcript> gene, double minReadsPerGene) {
		Gene plus = getGene(geneBase,useChrPrefix,firstID,gene,'+',minReadsPerGene);
		if(plus != null) {
			allGenes.add(plus);
			firstID++;
		}
		Gene minus = getGene(geneBase,useChrPrefix,firstID,gene,'-',minReadsPerGene);
		if(minus != null) {
			allGenes.add(minus);
			firstID++;
		}
		Gene none = getGene(geneBase,useChrPrefix,firstID,gene,'.',minReadsPerGene);
		if(none != null) {
			allGenes.add(none);
			firstID++;
		}
		return firstID;
	}

	private int len(int[] diff) {
		int l=0;
		while(l < diff.length && diff[l]> -1) {
			l++;
		}
		return l;
	}
	
	
	
	public LinkedList<Transcript> quantify(LinkedList<Transcript> original, double percentExplained, double minReadsPerTranscript, double maxFraction, 
			double percentAbundance, double scaleIntronReads, double delta, int nIterations, double scale,
			int minIntronLength, boolean longReads) throws CloneNotSupportedException{
		
		if(original.size() == 0) {
			return original;
		}
		
		boolean[][] inoutExons = new boolean[original.size()][contRegs.length];
		boolean[][] inoutIntrons = new boolean[original.size()][introns.length];
		
		int i=0;
		for(Transcript t : original) {
			for(Integer j : t.exons) {
				inoutExons[i][j] = true;
			}
			for(Integer j : t.introns) {
				inoutIntrons[i][j] = true;
			}
			i++;
		}
		
		int totalNumReads = 0;
		for(i=0;i<contRegs.length;i++) {
			int num = contRegs[i].readIdxs[ contRegs[i].readIdxs.length-1 ];
			if(num > totalNumReads) {
				totalNumReads = num;
			}
		}
		for(i=0;i<introns.length;i++) {
			if(introns[i].readIdxs.length==0) {
				System.out.println(Arrays.toString(contRegs));
				System.out.println(Arrays.toString(introns));
				System.out.println(introns[i].edge.getRelStart()+" "+introns[i].edge.getRelEnd()+" "+regionStart);
			}
			int num = introns[i].readIdxs[ introns[i].readIdxs.length-1 ];
			if(num > totalNumReads) {
				totalNumReads = num;
			}
		}
		
		if(totalNumReads == 0) {
			return new LinkedList<SplicingGraph.Transcript>();
		}else {
			totalNumReads++;
		}
		
		
		boolean[][] isAlternative = new boolean[contRegs.length][introns.length];
		for(i=0;i<contRegs.length;i++) {
			int startContReg = contRegs[i].regionStart;
			int endContReg = contRegs[i].regionEnd;
			for(int j=0;j<introns.length;j++) {
				int intronStart = introns[j].getStart();
				int intronEnd = introns[j].getEnd();
				if(intronStart+1<intronEnd) {
					if( (intronStart >= startContReg && intronStart < endContReg) ||
							(startContReg >= intronStart && startContReg < intronEnd) ) {
						isAlternative[i][j] = true;
					}
				}
			}
		}
				
		double[] readWeights = new double[totalNumReads];
		double[] readNs = new double[totalNumReads];
		Arrays.fill(readWeights, 1.0);
		for(i=0;i<introns.length;i++) {
			if(introns[i].getStart()+1<introns[i].getEnd()) {
				double maxNum = 1;
				for(int k=0;k<contRegs.length;k++) {
					if(isAlternative[k][i]) {
						double num = contRegs[k].readIdxs.length;
						if(num>maxNum) {
							maxNum = num;
						}
					}
				}
				double fac = 2.0*introns[i].readIdxs.length/(maxNum+introns[i].readIdxs.length);
				
				for(int j=0;j<introns[i].readIdxs.length;j++) {
					readWeights[introns[i].readIdxs[j]] *= fac;
					readNs[introns[i].readIdxs[j]]++;
				}
			}
		}
		

		for(i=0;i<contRegs.length;i++) { 
			double maxNum = 1;
			for(int k=0;k<introns.length;k++) {
				if(isAlternative[i][k]) {
					double num = introns[k].readIdxs.length;
					if(num > maxNum) {
						maxNum = num;
					}
				}
			}
			double fac = 2.0*contRegs[i].readIdxs.length/(maxNum+contRegs[i].readIdxs.length);

			for(int j=0;j<contRegs[i].readIdxs.length;j++) {
				
				readWeights[contRegs[i].readIdxs[j]] *= fac;
				readNs[contRegs[i].readIdxs[j]]++;
			}
		}
		

		for(i=0;i<readWeights.length;i++) {
				if(readNs[i] > 0) {
					readWeights[i] = Math.pow(readWeights[i], 1.0/readNs[i]);
				}else {
					readWeights[i] = 0.0;
				}
		}

		if(longReads) {//long reads may cover *many* introns
			
			int[] intronNs = new int[readNs.length];
			for(i=0;i<introns.length;i++) {
				if(introns[i].getStart()+minIntronLength<introns[i].getEnd()) {
					for(int j=0;j<introns[i].readIdxs.length;j++) {
						intronNs[introns[i].readIdxs[j]]++;
					}
				}
			}
			int maxIN = 0;
			for(i=0;i<intronNs.length;i++) {
				if(intronNs[i] > maxIN) {
					maxIN = intronNs[i];
				}
			}
			
			
			double[] facs = new double[maxIN+1];
			facs[0] = 1.0;
			if(maxIN>0) {
				facs[1] = scaleIntronReads;
			}
			double rat = 1.5;
			for(i=2;i<facs.length;i++) {
				facs[i] = facs[i-1]*( 1.0 + (scaleIntronReads-1.0)/rat );
				rat *= 1.5;
			}
			
			for(i=0;i<readNs.length;i++) {
				readWeights[i] *= facs[intronNs[i]];
			}
			
		}else {
			for(i=0;i<introns.length;i++) {
				if(introns[i].getStart()+minIntronLength<introns[i].getEnd()) {
					for(int j=0;j<introns[i].readIdxs.length;j++) {
						readWeights[introns[i].readIdxs[j]] *= scaleIntronReads;
					}
				}
			}
		}
		
		boolean[][] readsTimesTranscripts = new boolean[totalNumReads][original.size()];
		
		for(i=0;i<original.size();i++) {
			
			int nUnion = 1;
			int nDiff = 1;
			for(int j=0;j<inoutExons[i].length;j++) {
				if(inoutExons[i][j]) {
					nUnion += contRegs[j].readIdxs.length;
				}else {
					nDiff += contRegs[j].readIdxs.length;
				}
			}
			for(int j=0;j<inoutIntrons[i].length;j++) {
				if(inoutIntrons[i][j]) {
					nUnion += introns[j].readIdxs.length;
				}else {
					nDiff += introns[j].readIdxs.length;
				}
			}
			
			
			
			int[] tempUnion = new int[nUnion];
			int[] union = new int[nUnion];
			union[0] = tempUnion[0] = -1;
			
			int[] tempDiff = new int[nDiff];
			int[] diff = new int[nDiff];
			diff[0] = tempDiff[0] = -1;
			
			
			for(int j=0;j<inoutExons[i].length;j++) {
				if(inoutExons[i][j]) {
					unionSortedArrays(tempUnion, union, contRegs[j].readIdxs);
					int[] temp = tempUnion;
					tempUnion = union;
					union = temp;
				}else {
					unionSortedArrays(tempDiff, diff, contRegs[j].readIdxs);
					int[] temp = tempDiff;
					tempDiff = diff;
					diff = temp;
				}
				
			}
			for(int j=0;j<inoutIntrons[i].length;j++) {
				if(inoutIntrons[i][j]) {
					unionSortedArrays(tempUnion, union, introns[j].readIdxs);
					int[] temp = tempUnion;
					tempUnion = union;
					union = temp;
				}else {
					unionSortedArrays(tempDiff, diff, introns[j].readIdxs);
					int[] temp = tempDiff;
					tempDiff = diff;
					diff = temp;
				}
				
			}
			
			int[] compatibleReads = setdiff(union,diff);
			
			for(int j=0;j<compatibleReads.length;j++) {
				readsTimesTranscripts[compatibleReads[j]][i] = true;
			}
			
		}
		double[] relTranscriptLen = new double[inoutExons.length];
		Arrays.fill(relTranscriptLen, 1.0);
		double[] intronWeights = new double[inoutIntrons.length];
		Arrays.fill(intronWeights, 1.0);
		
		double[] cdsWeights = getCDSWeights(original);
		
		
		
		boolean[][] rTT_back = readsTimesTranscripts;//ArrayHandler.clone(readsTimesTranscripts);
		double[] rW_back = readWeights;//.clone();
		
		Pair<boolean[][],double[]> pair = makeUnique(readsTimesTranscripts, readWeights);
		readsTimesTranscripts = pair.getFirstElement();
		readWeights = pair.getSecondElement();
		
		double nReads = 0.0;
		double[] aPrioriTranscripts = new double[original.size()];
		double[][] gamma = new double[readsTimesTranscripts.length][original.size()];
		for(i=0;i<gamma.length;i++) {
			boolean found = false;
			for(int j=0;j<gamma[i].length;j++) {
				if(readsTimesTranscripts[i][j]) {
					gamma[i][j] = 1.0/relTranscriptLen[j]*intronWeights[j]*cdsWeights[j];
					found = true;
				}
			}
			if(found) {
				Normalisation.sumNormalisation(gamma[i]);
				nReads += readWeights[i];
			}
		}



		double oldD = Double.NEGATIVE_INFINITY;
		i=0;
		while(true) {
			
			aPriori(aPrioriTranscripts,gamma, readWeights);
			
			double d = gamma(gamma,aPrioriTranscripts,readsTimesTranscripts, readWeights,relTranscriptLen, intronWeights, cdsWeights);
		
			double diff = d - oldD;
			if(diff < delta || i>nIterations) {
				break;
			}

			oldD = d;
			i++;
		}

		readsTimesTranscripts = rTT_back;
		readWeights = rW_back;
		
		
		
		
		double[] idxs = new double[original.size()];
		for(i=0;i<idxs.length;i++) {
			idxs[i] = i;
		}

		
		double[] aPrioriTranscriptsOriginalOrder = aPrioriTranscripts.clone();
				
		ToolBox.sortAlongWith(aPrioriTranscripts, idxs);
		
		double sum = 0;
		i=aPrioriTranscripts.length;
		while(i>0 && aPrioriTranscripts[i-1]*nReads >= minReadsPerTranscript) {
			i--;
			sum += aPrioriTranscripts[i];
			if(sum >= percentExplained) {
				break;
			}
			if(i > 0) {
				
				if(aPrioriTranscripts[i]/aPrioriTranscripts[i-1] > maxFraction) {
					break;
				}
				if(aPrioriTranscripts[i-1]<percentAbundance) {
					break;
				}
			}
		}
		
		LinkedList<Transcript> result = new LinkedList<SplicingGraph.Transcript>();
		
		for(int j=idxs.length-1;j>=i;j--) {
			Transcript temp = original.get((int)idxs[j]);
			temp.setAbundance(aPrioriTranscripts[j]*nReads);
			temp.setStrand( getStrand( aPrioriTranscriptsOriginalOrder, readsTimesTranscripts, readWeights, (int)idxs[j], relTranscriptLen, intronWeights, cdsWeights ) );
			result.add(temp);
		}
		
		return result;
		
	}
	
	
	private double[] getCDSWeights(LinkedList<Transcript> original) {
		double[] w = new double[original.size()];
		
		double sum = 0.0;
		for(int i=0;i<w.length;i++) {
			int len = original.get(i).cdsEnd-original.get(i).cdsStart+1;
			w[i] = len;
			sum += len;
		}
		
		sum /= w.length;
		
		for(int i=0;i<w.length;i++) {
			w[i] /= sum;
		}
		
		return w;
		
	}

	private double[] getRelativeTranscriptLengths(boolean[][] inoutExons) {
		double[] rel = new double[inoutExons.length];
		for(int i=0;i<rel.length;i++) {
			double len = 0;
			for(int j=0;j<inoutExons[i].length;j++) {
				if(inoutExons[i][j]) {
					len += contRegs[j].regionEnd - contRegs[j].regionStart;
				}
			}
			rel[i] = len;
		}
		double mean = ToolBox.mean(rel);
		for(int i=0;i<rel.length;i++) {
			rel[i] = rel[i]/mean;//mean/rel[i];
		}
		Arrays.fill(rel, 1.0);
		return rel;
	}
	
	private double[] getNumberOfIntrons(boolean[][] inoutIntrons, int minIntronLength) {
		double[] num = new double[inoutIntrons.length];
		for(int i=0;i<num.length;i++) {
			int n = 0;
			for(int j=0;j<inoutIntrons[i].length;j++) {
				if(inoutIntrons[i][j] && introns[j].getStart() + minIntronLength < introns[j].getEnd()) {
					n++;
				}
			}
			num[i] = n;
		}
		return num;
	}
	
	private double[] getIntronWeights(boolean[][] inoutIntrons, int minIntronLength) {
		double[] w = getNumberOfIntrons(inoutIntrons, minIntronLength);
		double[] w2 = new double[w.length];
		for(int i=0;i<w2.length;i++) {
			w2[i] = getIntronWeight(w, i);
		}
		return w2;
	}

	private static Pair<boolean[][],double[]> makeUnique(boolean[][] readsTimesTranscripts, double[] readWeights2) throws CloneNotSupportedException {
		
		if(readsTimesTranscripts.length == 0) {
			return new Pair<boolean[][], double[]>(readsTimesTranscripts, readWeights2);
		}
		if(readsTimesTranscripts.length != readWeights2.length) {
			throw new RuntimeException();
		}
		
		Integer[] sortIdxs = new Integer[readsTimesTranscripts.length];
		for(int i=0;i<sortIdxs.length;i++) {
			sortIdxs[i] = i;
		}
		
		Arrays.sort(sortIdxs,new Comparator<Integer>(){
			@Override
			public int compare(Integer o1, Integer o2) {
				for(int i=0;i<readsTimesTranscripts[o1].length;i++) {
					if(readsTimesTranscripts[o1][i] && !readsTimesTranscripts[o2][i]) {
						return 1;
					}else if(!readsTimesTranscripts[o1][i] && readsTimesTranscripts[o2][i]) {
						return -1;
					}
				}
				return 0;
			}
		});		
		
		int uni = 1;
		for(int i=1;i<sortIdxs.length;i++) {
			boolean eq = true;
			boolean[] x = readsTimesTranscripts[sortIdxs[i]];
			boolean[] y = readsTimesTranscripts[sortIdxs[i-1]];
			for(int j=0;j<x.length;j++) {
				if(x[j] != y[j]) {
					eq = false;
					break;
				}
			}
			if(!eq) {
				uni ++;
			}
		}
		
		boolean[][] rTT = new boolean[uni][];
		double[] rW = new double[uni];
		
		rTT[0] = readsTimesTranscripts[sortIdxs[0]];
		rW[0] = readWeights2[sortIdxs[0]];
		int idx = 0;
		
		for(int i=1;i<sortIdxs.length;i++) {
			boolean eq = true;
			boolean[] x = readsTimesTranscripts[sortIdxs[i]];
			boolean[] y = readsTimesTranscripts[sortIdxs[i-1]];
			for(int j=0;j<x.length;j++) {
				if(x[j] != y[j]) {
					eq = false;
					break;
				}
			}
			if(eq) {
				rW[idx] += readWeights2[sortIdxs[i]];
			}else {
				idx++;
				rTT[idx] = x;
				rW[idx] = readWeights2[sortIdxs[i]];
			}
		}		
		
		return new Pair<boolean[][],double[]>(rTT,rW);
	}
	
	private void printPercentUnique(boolean[][] readsTimesTranscripts) {
		try {
			boolean[][] x = ArrayHandler.clone(readsTimesTranscripts);
			
			Arrays.sort(x,new Comparator<boolean[]>(){
				@Override
				public int compare(boolean[] o1, boolean[] o2) {
					for(int i=0;i<o1.length;i++) {
						if(o1[i] && !o2[i]) {
							return 1;
						}else if(!o1[i] && o2[i]) {
							return -1;
						}
					}
					return 0;
				}
			});
			
			int uni = 1;
			for(int i=1;i<x.length;i++) {
				boolean eq = true;
				for(int j=0;j<x[i].length;j++) {
					if(x[i][j] != x[i-1][j]) {
						eq = false;
						break;
					}
				}
				if(!eq) {
					uni ++;
				}
			}
			
			System.out.println("%: "+uni+" "+x.length);
			
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}		
	}

	private double delta(double[] aPrioriPrev, double[] aPrioriTranscripts) {
		double d = 0;
		for(int i=0;i<aPrioriPrev.length;i++) {
			d += (aPrioriPrev[i] - aPrioriTranscripts[i])*(aPrioriPrev[i] - aPrioriTranscripts[i]);
		}
		return d;
	}

	private char getStrand(double[] aPrioriTranscripts, boolean[][] readsTimesTranscripts, double[] readWeights, int t, double[] relTranscriptLen, double[] intronWeights, double[] cdsWeights) {
		double[] counts = new double[3];
		for(int i=0;i<readsTimesTranscripts.length;i++) {
			char strand = this.strands.get(i);
			
			double[] gamma = new double[aPrioriTranscripts.length];
			gamma(gamma,aPrioriTranscripts,readsTimesTranscripts[i],readWeights[i], relTranscriptLen, intronWeights, cdsWeights);
			
			if(strand == '+') {
				counts[0] += gamma[t];
			}else if(strand == '-') {
				counts[2] += gamma[t];
			}else {
				counts[1] += gamma[t];
			}
		}
		Normalisation.sumNormalisation(counts);
		
		if(counts[0] > 0.8) {
			return '+';
		}else if(counts[2] > 0.8) {
			return '-';
		}else {
			return '.';		
		}
	}
	
	private static final double gamma(double[] gamma, double[] aPrioriTranscripts, boolean[] readsTimesTranscripts, double readWeight, double[] relTranscriptLen, double[] intronWeights, double[] cdsWeights) {
		boolean found = false;
		for(int j=0;j<gamma.length;j++) {
			if(readsTimesTranscripts[j]) {
				gamma[j] = aPrioriTranscripts[j] / relTranscriptLen[j] * intronWeights[j] * cdsWeights[j];
				found = true;
			}
		}
		
		if(found) {
			double temp = sumNormalisation(gamma);
			return readWeight*Math.log( temp );
		}else {
			return 0;
		}
	}
	

	private static final double gamma(double[][] gamma, double[] aPrioriTranscripts, boolean[][] readsTimesTranscripts, double[] readWeights, double[] relTranscriptLen, double[] intronWeights, double[] cdsWeights) {
		double d = 0.0;
		for(int i=0;i<gamma.length;i++) {
			
			d += gamma(gamma[i],aPrioriTranscripts,readsTimesTranscripts[i], readWeights[i], relTranscriptLen, intronWeights, cdsWeights);
			
		}
		return d;
	}
	
	private static double getIntronWeight(double[] numIntrons, int j) {
		/*double ma = ToolBox.max(numIntrons);
		if(numIntrons[j] == ma) {
			return 1.5;
		}else {
			return 1.0 + 0.3*numIntrons[j]/ma;
		}*/
		return 1.0;
	}
	
	public static final double sumNormalisation(double[] gamma) {
		double sum = 0.0;
		for(int i=0;i<gamma.length;i++) {
			sum += gamma[i];
		}
		for(int i=0;i<gamma.length;i++) {
			gamma[i] /= sum;
		}
		return sum;
	}

	private static final void aPriori(double[] aPrioriTranscripts, double[][] gamma, double[] readWeights) {
		Arrays.fill(aPrioriTranscripts, 0.0);
		for(int i=0;i<gamma.length;i++) {
			for(int j=0;j<gamma[i].length;j++) {
				aPrioriTranscripts[j] += gamma[i][j]*readWeights[i];
			}
		}
		sumNormalisation(aPrioriTranscripts);
	}
	
	

	public LinkedList<Transcript> optimize(LinkedList<Transcript> original, double lambda1, double lambda2) throws DimensionException, TerminationException, IOException, EvaluationException, Exception {
		
		
		
		double[] target = new double[contRegs.length];
		for(int i=0;i<target.length;i++) {
			target[i] = (double)contRegs[i].readIdxs.length;
		}
		
		System.out.println(Arrays.toString(target));
		
		int[][] assignment = new int[original.size()][];
		int i=0;
		for(Transcript t : original) {
			assignment[i] = new int[t.exons.size()];
			int j=0;
			for(Integer k : t.exons) {
				assignment[i][j] = k;
				j++;
			}
			i++;
		}
		
		
		double[] w = new double[original.size()+1];
		Arrays.fill(w, 1.0);
		
		TranscriptLasso lasso = new TranscriptLasso(assignment, target, lambda1, lambda2);
		
		
		Optimizer.optimize(Optimizer.CONJUGATE_GRADIENTS_PRP, lasso, w, new SmallDifferenceOfFunctionEvaluationsCondition(1E-9), 1E-9, new ConstantStartDistance(1E-4), SafeOutputStream.getSafeOutputStream(System.out));
		
		System.out.println(Arrays.toString(w));
		
		LinkedList<Transcript> opt = new LinkedList<SplicingGraph.Transcript>();
		
		double thresh = Math.min(0.0, ToolBox.max(w));
		
		i=0;
		for(Transcript t : original) {
			if(w[i]>=thresh) {
				opt.add(t);
			}
			i++;
		}
		
		return opt;
		
	}
	
	public void enumerateTranscripts(LinkedList<Transcript> results, double minReads, double minFraction) {
		
		LinkedList<ComparableElement<Integer, Double>> starts = new LinkedList<ComparableElement<Integer,Double>>();
		for(int i=0;i<contRegs.length;i++) {
			if(contRegs[i].left.size() == 0) {
				starts.add( new ComparableElement<Integer, Double>(i, -(double)contRegs[i].readIdxs.length) );
			}
		}
		Collections.sort(starts);
		for(ComparableElement<Integer, Double> start : starts) {
		//	System.out.println(start);
			LinkedList<Integer> exons = new LinkedList<Integer>();
			exons.add(start.getElement());
			LinkedList<Integer> introns = new LinkedList<Integer>();
			enumerate( exons, introns, results, minReads, minFraction );
		}
		
	}
	
	
	private void enumerate(LinkedList<Integer> exons, LinkedList<Integer> introns, LinkedList<Transcript> results, double minReads, double minFraction) {
		if(results.size() > 1000) {
			System.out.println("#"+results.size());
			return;
		}
		if( adjList[exons.getLast()].size() == 0 ) {
			if(contRegs[ exons.getLast() ].readIdxs.length>=minReads) {
				Transcript temp = new Transcript(exons,introns);
				results.add(temp);
			}
		}else {
			boolean rec = false;
			int forLater = -1;
			int numForLater = 0;
			for(ComparableElement<Integer,Integer> liEl : adjList[exons.getLast()]) {
				if(this.introns[ liEl.getWeight() ].readIdxs.length >= minReads
						&& contRegs[ liEl.getElement() ].readIdxs.length >= minReads) {

						double max = Math.max( contRegs[ exons.getLast() ].readNodes.getLast().getNumberOfReads() , contRegs[liEl.getElement()].readNodes.getFirst().getNumberOfReads() );
						double inter = this.introns[liEl.getWeight() ].readIdxs.length;

						if(max*minFraction < inter) {

							exons.addLast(liEl.getElement());
							introns.addLast(liEl.getWeight());
													
							enumerate(exons, introns, results, minReads, minFraction);
							
							rec = true;
							
							exons.removeLast();
							introns.removeLast();

						}
				}else if(contRegs[ liEl.getElement() ].readIdxs.length >= minReads) {
					if(contRegs[ liEl.getElement() ].readIdxs.length > numForLater) {
						numForLater = contRegs[ liEl.getElement() ].readIdxs.length;
						forLater = liEl.getElement();
					}
				}
			}
			
		}
		
	}
	
	
	public void enumerateTranscripts2(LinkedList<Transcript> results, double minReads, double minFraction, double maxNum) {
		
		LinkedList<ComparableElement<Integer, Double>> starts = new LinkedList<ComparableElement<Integer,Double>>();
		

		double totalNum = countTranscripts(minReads, minFraction);
		
		double totalSum = 0.0;
		for(int i=0;i<contRegs.length;i++) {
			if(contRegs[i].left.size() == 0) {
				starts.add( new ComparableElement<Integer, Double>(i, -(double)contRegs[i].readIdxs.length) );
				totalSum += contRegs[i].getNumberOfReads()*ToolBox.sum(adjCounts[i]);
			}
		}
		
		
		
		
		if(totalNum < maxNum) {
			enumerateTranscripts(results, minReads, minFraction );
			return;
		}
		
		double theoCount = 0;
		
		Collections.sort(starts);
		for(ComparableElement<Integer, Double> start : starts) {
			
			double prob = (double)(ToolBox.sum(adjCounts[start.getElement()])*contRegs[start.getElement()].getNumberOfReads())/totalSum;
			theoCount += prob*maxNum;			
			
			LinkedList<Integer> exons = new LinkedList<Integer>();
			exons.add(start.getElement());
			LinkedList<Integer> introns = new LinkedList<Integer>();
			
			enumerate2( exons, introns, results, minReads, minFraction, maxNum*prob );
			if(results.size() > maxNum) {
				break;
			}
			if(theoCount > results.size()) {
				maxNum += theoCount - results.size();
			}
		}
		
		
	}
	
	
	private void enumerate2(LinkedList<Integer> exons, LinkedList<Integer> introns, LinkedList<Transcript> results, double minReads, double minFraction, double maxTranscripts) {

		if( adjList[exons.getLast()].size() == 0 ) {
			if(contRegs[ exons.getLast() ].readIdxs.length>=minReads) {
				Transcript temp = new Transcript(exons,introns);
				results.add(temp);
			}
		}else {
						
			double totalSum = 0.0;
			for(ComparableElement<Integer,Integer> liEl : adjList[exons.getLast()]) {
				totalSum += contRegs[liEl.getElement()].getNumberOfReads()*adjCounts[exons.getLast()][liEl.getElement()];
			}
			
			if(totalSum ==0) {
				return;
			}
			
			
			
			double initial = results.size();
			double theoCount = 0.0;

			boolean rec = false;
			int forLater = -1;
			int numForLater = 0;
			for(ComparableElement<Integer,Integer> liEl : adjList[exons.getLast()]) {
				if(this.introns[ liEl.getWeight() ].readIdxs.length >= minReads
						&& contRegs[ liEl.getElement() ].readIdxs.length >= minReads) {



					double max = Math.max( contRegs[ exons.getLast() ].readNodes.getLast().getNumberOfReads() , contRegs[liEl.getElement()].readNodes.getFirst().getNumberOfReads() );
					double inter = this.introns[liEl.getWeight() ].readIdxs.length;

					if(max*minFraction < inter) {

						double prob = (double)(adjCounts[exons.getLast()][liEl.getElement()]*contRegs[liEl.getElement()].getNumberOfReads())/totalSum;

						theoCount += maxTranscripts*prob;

						exons.addLast(liEl.getElement());
						introns.addLast(liEl.getWeight());

						rec = true; 
						enumerate2(exons, introns, results, minReads, minFraction, maxTranscripts*prob );

						exons.removeLast();
						introns.removeLast();

						if(results.size() - initial > maxTranscripts) {
							break;
						}

						if(theoCount > results.size() - initial) {
							maxTranscripts += theoCount - (results.size() - initial);
						}

					}


				}else if(contRegs[ liEl.getElement() ].readIdxs.length >= minReads) {
					if(contRegs[ liEl.getElement() ].readIdxs.length > numForLater) {
						numForLater = contRegs[ liEl.getElement() ].readIdxs.length;
						forLater = liEl.getElement();
					}
				}
			}
			
		}

	}

	
	
	
	
	
	
	
	
	
	public double countTranscripts(double minReads, double minFraction) {
		
		LinkedList<ComparableElement<Integer, Double>> starts = new LinkedList<ComparableElement<Integer,Double>>();
		for(int i=0;i<contRegs.length;i++) {
			if(contRegs[i].left.size() == 0) {
				starts.add( new ComparableElement<Integer, Double>(i, -(double)contRegs[i].readIdxs.length) );
			}
		}
		double sum = 0;
		Collections.sort(starts);
		for(ComparableElement<Integer, Double> start : starts) {
			LinkedList<Integer> exons = new LinkedList<Integer>();
			exons.add(start.getElement());
			LinkedList<Integer> introns = new LinkedList<Integer>();
			sum += count( exons, introns, minReads, minFraction );
			if(sum>1E6) {
				return sum;
			}
		}
		
		return sum;
	}
	
	private double count(LinkedList<Integer> exons, LinkedList<Integer> introns, double minReads, double minFraction) {
		double count = 0.0;
		if( exons.size() > 0 && adjList[exons.getLast()].size() == 0 ) {
			if(contRegs[ exons.getLast() ].readIdxs.length>=minReads) {
				count++;
			}
		}else {
			
			for(ComparableElement<Integer,Integer> liEl : adjList[exons.getLast()]) {
				if(this.introns[ liEl.getWeight() ].readIdxs.length >= minReads
						&& contRegs[ liEl.getElement() ].readIdxs.length >= minReads) {
					
					double max = Math.max( contRegs[ exons.getLast() ].readNodes.getLast().getNumberOfReads() , contRegs[liEl.getElement()].readNodes.getFirst().getNumberOfReads() );
					double inter = this.introns[liEl.getWeight() ].readIdxs.length;
					
					if(max*minFraction < inter) {

						int temp = exons.getLast();
						
						exons.addLast(liEl.getElement());
						introns.addLast(liEl.getWeight());
						
						double temp2 = count(exons, introns, minReads, minFraction);
						
						adjCounts[temp][liEl.getElement()] += temp2;
						count += temp2;
						if(count>1E6) {
							return count;
						}
						exons.removeLast();
						introns.removeLast();
					}
				}
			}
		}
		return count;
	}


	private ContiguousRegion[] getContRegs(LinkedList<Integer> exons) {
		ContiguousRegion[] res = new ContiguousRegion[exons.size()];
		int j=0;
		for(int i : exons) {
			res[j] = contRegs[i];
			j++;
		}
		return res;
	}

	
	public class Gene {
		
		private LinkedList<Transcript> transcripts;
		private String id;
		private char strand;
		
		public Gene(String id, LinkedList<Transcript> transcripts, char strand) {
			this.id = id;
			this.transcripts = transcripts;
			this.strand = strand;
		}
		
		public String toString() {
			StringBuffer buf = new StringBuffer();
			
			int start = Integer.MAX_VALUE;
			int end = -1;
			
			for(Transcript t : transcripts) {
				if(t.getStart()<start) {
					start = t.getStart();
				}
				if(t.getEnd()>end) {
					end = t.getEnd();
				}
			}
			
			buf.append(
					chromosome+"\t"
					+ "GeMoRNA\t"
					+ "gene\t"
					+ start+"\t"
					+ end+"\t"
					+ ".\t"
					+ strand+"\t"
					+ ".\t"
					+ "ID="+id+";"
					+ "\n"
					);
			
			for(Transcript t : transcripts) {
				buf.append(t.toString());
			}
					
			return buf.toString();
		}
		
		public String getChrom() {
			return chromosome;
		}
		
	}
	

	public class Transcript implements Comparable{
		
		private LinkedList<Integer> exons;
		private LinkedList<Integer> introns;
		private String id;
		private String geneID;
		private double abundance;
		private char strand;
		private int cdsStart;
		private int cdsEnd;
		
		private int extraLeft;
		private int extraRight;
		
		private HashMap<String,String> attributes;
		
		private Transcript(LinkedList<Integer> exons, LinkedList<Integer> introns) {
			this.exons = (LinkedList<Integer>) exons.clone();
			this.introns = (LinkedList<Integer>) introns.clone();
			this.id = "T";
			this.strand = '.';
			this.cdsStart = -1;
			this.cdsEnd = -1;
		}
		
		//not used
		public char determineStrandFromReads() {
			int countplus = 0;
			int countminus = 0;
			for(int i=0;i<exons.size();i++) {
				int[] idxs = contRegs[exons.get(i)].readIdxs;
				for(int j=0;j<idxs.length;j++) {
					char temp = strands.get(j);
					if(temp == '+') {
						countplus++;
					}else if(temp == '-'){
						countminus++;
					}
				}
			}
			for(int i=0;i<introns.size();i++) {
				int[] idxs = SplicingGraph.this.introns[introns.get(i)].readIdxs;
				for(int j=0;j<idxs.length;j++) {
					char temp = strands.get(j);
					if(temp == '+') {
						countplus++;
					}else if(temp == '-'){
						countminus++;
					}
				}
			}
			if(countplus > countminus*2) {
				return '+';
			}else if(countminus > countplus*2) {
				return '-';
			}else {
				return '.';
			}
		}
		
		
		public double getAbundance() {
			return abundance;
		}

		public char getStrand() {
			return this.strand;
		}
		
		public void setStrand(char strand) {
			this.strand = strand;
		}
		
		public void setGene(String geneID) {
			this.geneID = geneID;
		}

		public String getId() {
			return this.id;
		}
		
		public void setId(String id) {
			this.id = id;
		}
		
		public void setAbundance(double abundance) {
			this.abundance = abundance;
		}
		
		public void setCDSStartEnd(int start, int end) {
			int l = getLength();
			while(start<0) {
				start +=3;
			}
			while(end>l) {
				end -= 3;
			}
			this.cdsStart = start;
			this.cdsEnd = end;
		}
		
		public String toString() {
			StringBuffer buf = new StringBuffer();
			
						
			int ce = 0;
			
			int i=0;
			int start = -1;
			int end = -1;
			int currPos = 0;
			//int phase = -1;
			
			char startAS = 'X';
			char stopAS = 'X';

			String[] exStrings = new String[exons.size()];
			String[][] cdsStrings = new String[exons.size()][];
			
			int[] lens = new int[exons.size()];
			
			while(i<exons.size()) {
				
				start = contRegs[exons.get(i)].regionStart;
				if(i==0) {
					start += extraLeft;
				}
				end = contRegs[exons.get(i)].regionEnd;
				
				while( i+1<exons.size() && contRegs[exons.get(i)].regionEnd+1 == contRegs[exons.get(i+1)].regionStart ) {
					i++;
					end = contRegs[exons.get(i)].regionEnd;
				}
				
				if(i==exons.size()-1) {
					end += extraRight;
				}
				
				int currLen = end-start+1;
				
				exStrings[i] = 
						chromosome+"\t"
						+ "GeMoRNA\t"
						+ "exon\t"
						+ (regionStart + start)+"\t"
						+ (regionStart + end)+"\t"
						+ ".\t"
						+ this.strand+"\t"
						+ ".\t"
						+"Parent="+id
						+ "\n";
				
				if(cdsStart > -1 && currPos + currLen > cdsStart && currPos <= cdsEnd) {
					
					int locCdsStart = currPos > cdsStart ? start : start + (cdsStart-currPos);
					int locCdsEnd = currPos + currLen <= cdsEnd ? end : start + (cdsEnd-(currPos));
					
					
					
					
					
					cdsStrings[i] = new String[] {
							chromosome+"\t"
							+ "GeMoRNA\t"
							+ "CDS\t"
							+ (regionStart + locCdsStart)+"\t"
							+ (regionStart + locCdsEnd)+"\t"
							+ ".\t"
							+ this.strand+"\t",
							"\t"
							+"Parent="+id
							+ "\n"
					};
					ce++;
					
					lens[i] = (locCdsEnd-locCdsStart+1);
					
				}
				
				
				currPos += currLen;
				
				
				i++;
			}
			
			
			double tie = 0.0;
			int minSplitReads = Integer.MAX_VALUE;
			if(SplicingGraph.this.introns.length==0 || SplicingGraph.this.introns[0].readIdxs == null) {
				tie = Double.NaN;
				minSplitReads = 0;
			}else {
				for(i=0;i<introns.size();i++) {
					int temp = SplicingGraph.this.introns[introns.get(i)].readIdxs.length;
					tie +=  temp > 0 ? 1 : 0;
					if(temp < minSplitReads) {
						minSplitReads = temp;
					}
				}
				tie /= introns.size();
			}
			
			int aa = 0;
			if(cdsStart > -1) {
				aa = getCDSLength()/3;
			}
			
			
			
			if(cdsStart>-1) {
				if(this.strand=='-') {
					String key = getCodon(cdsStart);
					
					key = Tools.rc(key);
					Character asChar = Genome.code.get(key);
					if(asChar != null) {
						stopAS = asChar;
					}
					key = getCodon(cdsEnd-2);
					
					key = Tools.rc(key);
					
					asChar = Genome.code.get(key);
					if(asChar != null) {
						startAS = asChar;
					}
				}else {
					String key = getCodon(cdsStart);
					
					Character asChar = Genome.code.get(key);
					if(asChar != null) {
						startAS = asChar;
					}
					key = getCodon(cdsEnd-2);
					
					asChar = Genome.code.get(key);
					if(asChar != null) {
						stopAS = asChar;
					}
				}				
			
			}
			
			int minCov = Integer.MAX_VALUE;
			double avgCov = 0;
			double n = 0;
			double tpc = 0;
			for(i=0;i<exons.size();i++) {
				for(Node node : contRegs[exons.get(i)].readNodes) {
					int temp = node.getNumberOfReads();
					if(temp < minCov) {
						minCov = temp;
					}
					avgCov += temp;
					if(temp > 0) {
						tpc++;
					}
					n++;
				}
			}
			for(i=0;i<introns.size();i++) {
				int temp = 0;
				if(SplicingGraph.this.introns.length!=0 && SplicingGraph.this.introns[0].readIdxs != null) {
					temp = SplicingGraph.this.introns[introns.get(i)].readIdxs.length;
				}
				if(temp < minCov) {
					minCov = temp;
				}
				avgCov += temp;
				if(temp > 0) {
					tpc++;
				}
				n++;
			}
			
			avgCov /= n;
			tpc /=n;
			
			
			buf.append(
					chromosome+"\t"
					+ "GeMoRNA\t"
					+ "mRNA\t"
					+ (regionStart + contRegs[exons.getFirst()].regionStart + extraLeft)+"\t"
					+ (regionStart + contRegs[exons.getLast()].regionEnd + extraRight)+"\t"
					+ ".\t"
					+ this.strand+"\t"
					+ ".\t"
					+ "ID="+id+";"
					+ "Parent="+geneID+";"
					//+ "abundance="+abundance+"; "
					+ "percentCanonical="+percentCanonicalSpliceSites()+";"
					+ "aa="+aa+";"
					+ "score="+df.format(abundance)+";"
					+ "ce="+ce+";"
					+ "minCov="+minCov+";"
					+ "avgCov="+df.format(avgCov)+";"
					+ "tpc="+tpc+";"
					+ "tie="+tie+";"
					+ "start="+startAS+";"
					+ "stop="+stopAS+";"
					+ "minSplitReads="+minSplitReads+";"
					+ (attributes == null ? "" : attributes.entrySet().stream().map( e -> e.getKey()+"="+e.getValue()).reduce("", (t,u) -> t+"; "+u)
							)
					+ "\n"
					);
			
			if(strand == '-') {
				int phase = 0;
				for(i=exStrings.length-1;i>=0;i--) {
					if(exStrings[i] != null) {
						buf.append(exStrings[i]);
					}
					if(cdsStrings[i] != null) {
						buf.append(cdsStrings[i][0]+phase+cdsStrings[i][1]);
						int rest = (3-phase + lens[i])%3;
						phase = (3-rest)%3;
					}
				}
			}else {
				int phase = 0;
				for(i=0;i<exStrings.length;i++) {
					if(exStrings[i] != null) {
						buf.append(exStrings[i]);
					}
					if(cdsStrings[i] != null) {
						buf.append(cdsStrings[i][0]+phase+cdsStrings[i][1]);
						int rest = (3-phase + lens[i])%3;
						phase = (3-rest)%3;
					}
				}
			}
			
			return buf.toString();
		}
		
		
		private String getCodon(int codonStart) {
			char[] chromosome = Genome.genome.getChromosome(SplicingGraph.this.chromosome);
			
			char[] codon = new char[3];
			int ci = 0;
			
			int i=0;
			int off = 0;
			for(Integer idx : exons) {
				
				int start = regionStart + contRegs[idx].regionStart - 1;
				int end = regionStart + contRegs[idx].regionEnd;
				if(i==0) {
					start = regionStart + contRegs[idx].regionStart - 1 + extraLeft;
				}
				if(i==exons.size()-1) {
					end = regionStart + contRegs[idx].regionEnd + extraRight;
				}
				
				if(start<0) {
					start=0;
				}
				if(end>chromosome.length) {
					end = chromosome.length;
				}
				if(end>start) {
					
					if(off <= codonStart+ci && off+(end-start)>codonStart+ci) {
						for(int k=ci,l=0;k<3 && start+codonStart-off+k<end;k++,l++) {
							codon[k] = chromosome[start+codonStart-off+k];
							ci++;
						}
					}
					
					off += end-start;
				}
				i++;
			}
			return new String(codon);
		}
		
		
		private char[] getSequence(int extraLeft, int extraRight) {
			
			char[] chromosome = Genome.genome.getChromosome(SplicingGraph.this.chromosome);
			
			int l = this.getLength(extraLeft, extraRight);
			char[] res = new char[l];
			
			
			int i=0;
			int off = 0;
			for(Integer idx : exons) {
				
				int start = regionStart + contRegs[idx].regionStart - 1;
				int end = regionStart + contRegs[idx].regionEnd;
				if(i==0) {
					start = regionStart + contRegs[idx].regionStart - 1 + extraLeft;
				}
				if(i==exons.size()-1) {
					end = regionStart + contRegs[idx].regionEnd + extraRight;
				}
				
				if(start<0) {
					start=0;
				}
				if(end>chromosome.length) {
					end = chromosome.length;
				}
				if(end>start) {
					System.arraycopy(chromosome, start, res, off, end-start);
					
					off += end-start;
				}
				i++;
			}
			
			if(off < res.length) {
				char[] res2 = new char[off];
				System.arraycopy(res, 0, res2, 0, off);
				res = res2;
			}
			
			
			return res;
		}
		
		
		public int getLength() {
			return getLength(extraLeft,extraRight);
		}
		
		public int getLength(int extraLeft, int extraRight) {
			
			int l = 0;			
			for(int i=0;i<exons.size();i++) {
				int idx = exons.get(i);
				
				int start = contRegs[idx].regionStart - 1;
				int end = contRegs[idx].regionEnd;
				if(i==0) {
					start = contRegs[idx].regionStart - 1 + extraLeft;
				}
				if(i==exons.size()-1) {
					end = contRegs[idx].regionEnd + extraRight;
				}
				
				l += end-start;
				
				
				
			}
			
			return l;
		}
		
		
		public void addStrandAndCDS(int minProteinLength, boolean forTestSplit) {
			addStrandAndCDS(false,minProteinLength, forTestSplit);
		}
		
		public void addStrandAndCDS(boolean fix, int minProteinLength, boolean forTestSplit) {
						
			if(this.strand == '.') {
				this.strand = getStrandByIntrons();
			}
			
			int tl = this.getLength();
			if(tl + 2*Math.min(50, tl) < minProteinLength*3) {
				return;
			}
			
			char tempStrand = this.strand;
			
			
			if(fix) {
				this.extraLeft = this.extraRight = 0;
			}
			
			char[] seqOrig = this.getSequence(this.extraLeft,this.extraRight);
			char[] seqRc = null;
			char[] longSeqOrig = null;
			char[] longSeqRc = null;
			
			
			int maxLonger = Math.min(50, seqOrig.length/3);
			
			char[] chromosome = Genome.genome.getChromosome(SplicingGraph.this.chromosome);

			int tempLeft = -this.extraLeft;
			int tempRight = this.extraRight;
			if(fix) {
				tempLeft = Math.min((regionStart + contRegs[exons.getFirst()].regionStart - 1)/3, maxLonger)*3;
				tempRight = Math.min( (  chromosome.length - ( regionStart + contRegs[exons.getLast()].regionEnd) )/3, maxLonger)*3;
				longSeqOrig = this.getSequence(-tempLeft,tempRight);
			}
			
			HashMap<String,Character> code = Genome.code;

			
			int firstD = 1;
			int lastD = -1;
			if(this.strand == '+') {
				lastD = 1;
			}else if(this.strand == '-') {
				firstD = -1;
			}
			
			
			
			int numIntrons = getNumberOfIntrons();
			
			CDSCandidateList candList = new CDSCandidateList(3);
			StopIterator it = new StopIterator();
			
			for(int d=firstD;d>=lastD;d-=2) {
				
				char[] seq = seqOrig;
				char[] longSeq = longSeqOrig;
				
				if(d < 0) {
					seqRc = Tools.rc(seqOrig);
					seq = seqRc;
					if(longSeqOrig != null) {
						longSeqRc = Tools.rc(longSeqOrig);
						longSeq = longSeqRc;
					}
				}
				
				for(int i=1;i<=3;i++) {
					it.set(seq,i-1);
					
					boolean beforeContainsStart = false;
					boolean afterContainsStop = false;
					if(longSeq != null) {
						
						int lastStop = containsLastStop(longSeq, i-1, 0, (d<0 ? tempRight : tempLeft) + i-1 - 3) - (i-1);
						if(lastStop < 0) {
							lastStop = 0;
						}
						beforeContainsStart = containsStart(longSeq, i-1, lastStop, (d<0 ? tempRight : tempLeft) + (i-1)) >= 0;
						
						afterContainsStop = containsStop(longSeq,i-1,seq.length + (d<0 ? tempRight : tempLeft) - (seq.length-(i-1))%3 - (i-1)) >= 0;
						
					}					
					
					int off = 0;
					
					
					int j=0;
					while(it.hasNext()) {
						int[] lenstart = it.next();
						
						int length = lenstart[0];
						int start = lenstart[1];
						if(start < 0) {
							if(j==0 && beforeContainsStart) {
								start = 0;
							}else {
								start = length+1;
							}
						}
						int len = length - start + 1;
						if(it.hasNext() || it.isAtEnd() || afterContainsStop) {
							int numCDSIntrons = 0;
							if(d>0) {
								numCDSIntrons = getNumberOfCDSIntrons((off+start)*3,(off+length+1)*3);
							}else {
								numCDSIntrons = getNumberOfCDSIntrons(seq.length-(off+length+1)*3,seq.length-(off+start)*3);
							}
							candList.add(len, d*i, numCDSIntrons, numIntrons, j, !it.hasNext());
						}
						off += length+1;
						
						j++;
					}
				}
			}
						
			int[] bestCand = candList.getBest();
			
			int maxLen = bestCand[0];
			int frame = bestCand[1];
			int maxIdx = bestCand[2];
			boolean isLast = (bestCand[3] == 1);
			
			
			if(maxLen + tempLeft/3 + tempRight/3 < minProteinLength ) {
				return;
			}
			
			
			
			char[] seq = seqOrig;//this.getSequence(this.extraLeft,this.extraRight);
			char[] longSeq = longSeqOrig;
			if(frame < 0) {
				seq = seqRc;//Tools.rc(seq);
				longSeq = longSeqRc;
			}
			
			if(this.strand == '.') {
				this.strand = '+';
				if(frame < 0) {
					this.strand = '-';
				}
			}
			
			if(this.strand == '-') {
				frame = -frame;
			}
			if(frame != 0) {
				
				//String prot = Tools.translate(frame-1, new String(seq), code, false, Ambiguity.AMBIGUOUS);
				int protlen = (seq.length - (frame-1))/3;
				//System.out.println(prot);
				it.set(seq, frame-1);
				
			//	String[] parts = prot.split("\\*",-1);
				
				if(fix) {
					boolean fixed = false;
					if(isLast) {
						
						int posStart = (strand == '-' ? tempRight : tempLeft) - (seq.length-(frame-1))%3 - (frame-1);
						int cpos2 = containsStop(longSeq, frame-1, seq.length+posStart ) - (strand == '-' ? tempRight : tempLeft) - (frame-1);
						
						int cpos = cpos2/3;
						
						if(cpos >= 0) {
							cpos =  cpos*3 - protlen*3 + 3;
							if(strand == '+') {
								this.extraRight = cpos;
							}else if(strand == '-') {
								this.extraLeft = -cpos;
							}
							fixed = true;
						}
					}
					//boolean b1 = parts[0].indexOf("M")<0;
					int stopPos = containsStop(seq,frame-1,0);
					boolean b2 = containsStart(seq,frame-1,0,stopPos<0 ? seq.length : stopPos)<0;
					
					if(maxIdx == 0 && b2/*parts[0].indexOf("M")<0*/) {
						
						int lastStop = containsLastStop(longSeq, frame-1, 0, (strand == '-' ? tempRight : tempLeft) + frame-1 - 3) - (frame-1);
						int cpos = lastStop/3;
						
						if(cpos < 0) {
							cpos = 0;
						}
						if(lastStop < 0) {
							lastStop = 0;
						}
						int firstStart = containsStart(longSeq,frame-1,lastStop,(strand == '-' ? tempRight : tempLeft) + frame-1) - (frame-1);
						cpos = firstStart/3;
						
						if(cpos >= 0) {
							
							firstStart = (strand == '-' ? tempRight : tempLeft) + (frame==1 ? 3 : 0) - firstStart;
							cpos = firstStart;
							
							if(strand == '+') {
								this.extraLeft = -cpos;
							}else if(strand == '-') {
								this.extraRight = cpos;
							}
							fixed = true;
						}
						
					}
					
					if(fixed) {
						addStrandAndCDS(false,minProteinLength,forTestSplit);
						return;
					}
				}
				
				it.set(seq, frame-1);
				int j=0;
				int off = frame-1;
				while(it.hasNext() && j < maxIdx) {
					int[] lenstart = it.next();
					//System.out.println(prot+"\n"+Arrays.toString(lenstart));
					int length = lenstart[0];
					off += (length+1)*3;
					
					j++;
				}
				int[] lenstart = it.next();
				
				int temp2 = lenstart[1]*3;
				if(temp2 < 0) {
					if(maxIdx > 0) {
						throw new RuntimeException();
					}else {
						temp2 = 0;
					}
				}
				off += temp2;
				
				
				//System.out.println((maxLen*3)+" "+seq.length()+" "+cdsStart+" "+cdsEnd+" "+strand);
				if( (!forTestSplit && maxLen*3 < seq.length/20.0) || maxLen<minProteinLength) {
					this.cdsStart = this.cdsEnd = -1;
					this.strand = tempStrand;
				}else {
					if(this.strand == '+') {
						this.cdsStart = off;
						this.cdsEnd = off + maxLen*3 - 1;
					}else {
						this.cdsStart = seq.length - (off + maxLen*3);
						this.cdsEnd = seq.length - off -1;
					}
				}
			}
			
			
		}
		
		
		private int containsLastStop(char[] seq, int offset, int start, int end) {
			for(int j=end;j-3>=offset + start;j-=3) {
				if(Genome.isStop[seq[j]][seq[j+1]][seq[j+2]]) {
					return j;
				}
			}
			return -3;
		}
		
		private int containsStop(char[] seq, int offset, int start) {
			for(int j=offset + start;j+3<=seq.length;j+=3) {
				if(Genome.isStop[seq[j]][seq[j+1]][seq[j+2]]) {
					return j;
				}
			}
			return -3;
		}

		private int containsStart(char[] seq, int offset, int start, int end) {
			for(int j=offset + start;j+3<=end;j+=3) {
				if(Genome.isStart[seq[j]][seq[j+1]][seq[j+2]]) {
					return j;
				}
			}
			return -3;
		}
		
		
		public char getStrandByIntrons() {
			
			
			double[][] intronVotes = getIntronVotes();
			double nplus = ToolBox.sum(intronVotes[0]);
			double nminus = ToolBox.sum(intronVotes[1]);
			
			char strand = this.strand;
			
			if( (nplus+1)/(nplus+nminus+2) >=0.75 ) {
				strand = '+';
			}else if((nminus+1)/(nplus+nminus+2) >= 0.75) {
				strand = '-';
			}
			
			return strand;		
			
		}

		public int getStart() {
			return regionStart + SplicingGraph.this.contRegs[exons.getFirst()].regionStart + extraLeft;
		}
		
		public int getEnd() {
			return regionStart + SplicingGraph.this.contRegs[exons.getLast()].regionEnd + extraRight;
		}
		

		@Override
		public int compareTo(Object o) {
			if(o instanceof Transcript) {
				Transcript t = (Transcript)o;
				int cmpStart = Integer.compare(this.getStart(),t.getStart());
				
				if(cmpStart == 0) {
					return Integer.compare(this.getEnd(),t.getEnd());
				}else {
					return cmpStart;
				}
				
			}else {
				return -1;
			}
		}

		public int getNumberOfIntrons() {
			int n=0;
			Iterator<Integer> it = exons.iterator();
			Integer idx = it.next();
			while(it.hasNext()) {
				int lastEnd = contRegs[idx].regionEnd;
				idx = it.next();
				int start = contRegs[idx].regionStart;

				if(lastEnd+1<start) {
					n++;
				}
			}
			return n;
		}
		
		
		public int getNumberOfCDSIntrons(int cdsStart, int cdsEnd) {
			int n=0;
			Iterator<Integer> it = exons.iterator();
			Integer idx = it.next();
			int len = contRegs[idx].regionEnd - contRegs[idx].regionStart + 1;
			while(it.hasNext()) {
				int lastEnd = contRegs[idx].regionEnd;
				idx = it.next();
				int start = contRegs[idx].regionStart;

				if(lastEnd+1<start && len>cdsStart && len<=cdsEnd) {
					n++;
				}
				len += contRegs[idx].regionEnd - contRegs[idx].regionStart + 1;
			}
			return n;
		}
		
		
		
		public double percentCanonicalSpliceSites() {
			
			double[][] intronVotes = getIntronVotes();
			
			Iterator<Integer> it = exons.iterator();
			Integer idx = it.next();
			
			double n = 0;
			double canonical = 0;
			
			int i=1;
			while(it.hasNext()) {
				
				int lastEnd = contRegs[idx].regionEnd;
				idx = it.next();
				int start = contRegs[idx].regionStart;

				if(lastEnd+1<start) {
					
					if(this.strand == '-') {
						canonical += intronVotes[1][i];
					}else if(this.strand == '+') {
						canonical += intronVotes[0][i];
					}else {
						canonical += intronVotes[0][i] + intronVotes[1][i];
					}
					n += 2;
				}
				i++;
			}

			if(n==0) {
				return -1;
			}else {
				return canonical/n;
			}
			
		}
		
		
		private double[][] getIntronVotes() {
			char[] chromosome = Genome.genome.getChromosome(SplicingGraph.this.chromosome);
			
			
			double[][] intronVotes = new double[2][exons.size()];
			
			Iterator<Integer> it = exons.iterator();
			Integer idx = it.next();
			
			int i=1;

			while(it.hasNext()) {
				
				int lastEnd = contRegs[idx].regionEnd;
				idx = it.next();
				int start = contRegs[idx].regionStart;

				if(lastEnd+1<start) {

					//String left = chromosome.substring(regionStart + lastEnd, regionStart + lastEnd+2);
					char left1 = chromosome[regionStart+lastEnd];
					char left2 = chromosome[regionStart+lastEnd+1];
					//String right = chromosome.substring(regionStart + start - 3, regionStart + start-1);
					char right1 = chromosome[regionStart + start - 3];
					char right2 = chromosome[regionStart + start - 2];
					
					if(/*"GT".equals(left)*/
							'G' == left1 && 'T' == left2) {
						intronVotes[0][i] += 1;
					}else if(/*"GC".equals(left)*/
							'G' == left1 && 'C' == left2) {
						intronVotes[0][i] += 0.5;
					}else if(/*"CT".equals(left)*/
							'C' == left1 && 'T' == left2) {
						intronVotes[1][i] += 1;
					}
					
					if(/*"AG".equals(right)*/
							'A' == right1 && 'G' == right2) {
						intronVotes[0][i] += 1;
					}else if(/*"AC".equals(right)*/
							'A' == right1 && 'C' == right2) {
						intronVotes[1][i] += 1;
					}else if(/*"GC".equals(right)*/
							'G' == right1 && 'C' == right2) {
						intronVotes[1][i] += 0.5;
					}
				}
				i++;
			}
			
			return intronVotes;	
		}
		
		
		
		
		private int[][] getExonsIntronsAndLengths(int[] possibleSplitPositions) {
			int[] numExon = new int[exons.size()];
			int num = 0;
			for(int i=1;i<exons.size();i++) {
				if(contRegs[exons.get(i-1)].regionEnd+1 != contRegs[exons.get(i)].regionStart) {
					num++;
				}
				numExon[i] = num;
			}
			
			int[] possibleSplitExons = new int[possibleSplitPositions.length];
			int[] possibleSplitIntrons = new int[possibleSplitPositions.length];
			Arrays.fill(possibleSplitExons, -1);
			Arrays.fill(possibleSplitIntrons, -1);
			
			for(int i=0;i<possibleSplitPositions.length;i++) {
				int splitpos = possibleSplitPositions[i];
				for(int j=0;j<exons.size();j++) {
					if(contRegs[exons.get(j)].regionStart+regionStart < splitpos && contRegs[exons.get(j)].regionEnd+regionStart > splitpos) {
						possibleSplitExons[i] = j;
						break;
					}
				}
				for(int j=0;j<exons.size()-1;j++) {
					if(contRegs[exons.get(j)].regionEnd+regionStart < splitpos && contRegs[exons.get(j+1)].regionStart+regionStart >= splitpos) {
						possibleSplitIntrons[i] = j;
						break;
					}
				}
			}
			
			int[] cumSum = new int[exons.size()];
			for(int i=0;i<exons.size();i++) {
				if(i==0) {
					cumSum[i] = contRegs[exons.get(i)].regionEnd - contRegs[exons.get(i)].regionStart+1; 
				}else {
					cumSum[i] = cumSum[i-1] + contRegs[exons.get(i)].regionEnd - contRegs[exons.get(i)].regionStart+1;
				}
			}
			return new int[][] {possibleSplitExons,possibleSplitIntrons,cumSum};
		}
		
		
		
		
		
		
		
		
		private int[][] getExonsAndLengths(int[] possibleSplitPositions) {
			int[] numExon = new int[exons.size()];
			int num = 0;
			for(int i=1;i<exons.size();i++) {
				if(contRegs[exons.get(i-1)].regionEnd+1 != contRegs[exons.get(i)].regionStart) {
					num++;
				}
				numExon[i] = num;
			}
			
			int[] possibleSplitExons = new int[possibleSplitPositions.length];
			Arrays.fill(possibleSplitExons, -1);
			for(int i=0;i<possibleSplitPositions.length;i++) {
				int splitpos = possibleSplitPositions[i];
				for(int j=0;j<exons.size();j++) {
					if(contRegs[exons.get(j)].regionStart+regionStart < splitpos && contRegs[exons.get(j)].regionEnd+regionStart > splitpos) {
						if(numExon[ j ] > 0 && numExon[j]<numExon[numExon.length-1]) {
							possibleSplitExons[i] = j;
						}
						break;
					}
				}
			}
			
			int[] cumSum = new int[exons.size()];
			for(int i=0;i<exons.size();i++) {
				if(i==0) {
					cumSum[i] = contRegs[exons.get(i)].regionEnd - contRegs[exons.get(i)].regionStart+1; 
				}else {
					cumSum[i] = cumSum[i-1] + contRegs[exons.get(i)].regionEnd - contRegs[exons.get(i)].regionStart+1;
				}
			}
			return new int[][] {possibleSplitExons,cumSum};
		}
		
		public void testSplitByORF(LinkedList<Transcript> res, int[] possibleSplitPositions, boolean checkAllHaveCDS, int minProteinLength) {
			possibleSplitPositions = possibleSplitPositions.clone();
			
			LinkedList<Transcript> curr = new LinkedList<SplicingGraph.Transcript>();
			curr.add(this);
			
			outerloop:
			while(curr.size() > 0) {
				
				Transcript t = curr.removeFirst();
				
				char tempStrand = t.strand;
				
				t.addStrandAndCDS(minProteinLength,true);
				
				int[][] temp = t.getExonsIntronsAndLengths(possibleSplitPositions);
				int[] possibleSplitExons = temp[0];
				int[] possibleSplitIntrons = temp[1];
				int[] cumSum = temp[2];
				
				int idx = -1;
				boolean intron = false;
				
				for(int i=0;i<possibleSplitExons.length;i++) {
					if(possibleSplitIntrons[i] >= 0) {
						if( cumSum[possibleSplitIntrons[i]] < t.cdsStart ) {
							
							idx = i;
							intron = true;
							
						}else if(cumSum[ possibleSplitIntrons[i] ] > t.cdsEnd) {
							
							idx = i;
							intron = true;
							break;
							
						}
						
					}else if(possibleSplitExons[i] > 0 && possibleSplitExons[i] < t.exons.size()-1) {		
						if( ( cumSum[ possibleSplitExons[i]-1 ] < t.cdsStart) ) {
							
							idx = i;
							intron = false;
						}else if( cumSum[ possibleSplitExons[i] ] > t.cdsEnd ) {
							
							idx = i;
							intron = false;
							break;
							
						}
					}
				}
				
				if(idx > -1) {
					Transcript[] s = null;
					
					if(intron) {
						s = t.splitIntoTwo(possibleSplitIntrons[idx],true);
						
					}else {
						s = t.splitIntoTwo(possibleSplitExons[idx],false);

					}
					
					//TODO now here (see few lines below)
					possibleSplitPositions[idx] = -1;
					
					//TODO added splice criterion
					if(t.getNumberOfIntrons() > 3 && (s[0].getNumberOfIntrons()==0 || s[1].getNumberOfIntrons()==0) ) {
						curr.add(t);
						continue;
					}
					//TODO end added
					
					s[0].addStrandAndCDS(minProteinLength,true);
					s[1].addStrandAndCDS(minProteinLength,true);
					
					
					if(s[0].cdsStart == -1) {
						s[0].extraLeft = 0;
						s[0].addStrandAndCDS(true,minProteinLength,false);
					}
					
					if(s[1].cdsStart == -1) {
						s[1].extraRight = 0;
						s[1].addStrandAndCDS(true,minProteinLength,false);
					}
					
					if(!checkAllHaveCDS || (s[0].cdsStart>-1 && s[1].cdsStart>-1) ) {
						
						s[0].strand = '.';
						s[1].strand = '.';
						s[0].cdsStart = s[0].cdsEnd = -1;
						s[1].cdsStart = s[1].cdsEnd = -1;
						s[0].extraLeft = t.extraLeft;
						s[1].extraRight = t.extraRight;
						
						curr.add(s[0]);
						curr.add(s[1]);
						continue;
					}else {
						curr.add(t);
						continue;
					}
				}
				
				t.strand = tempStrand;
				t.cdsEnd = t.cdsStart = -1;
				
				res.add(t);
				
			}
						
		}
		
		
		
		public void testSplitByORF2(LinkedList<Transcript> res, int minProteinLength) {
			
			char oldstrand = this.strand;
			
			this.addStrandAndCDS(minProteinLength,true);
			
			int cdsStartIdx = -1;
			int cdsEndIdx = -1;
			
			int currPos = 0;
			
			int i=0;
			while(i<exons.size()) {
				
				int start = contRegs[exons.get(i)].regionStart;
				int end = contRegs[exons.get(i)].regionEnd;
				
				
				int currLen = end-start+1;
				
				if(cdsStart >= currPos && cdsStart < currPos+currLen ) {
					cdsStartIdx = i;
				}
				if(cdsEnd >= currPos && cdsEnd < currPos+currLen) {
					cdsEndIdx = i;
				}
				currPos += currLen;
				i++;
			}
			
		//	System.out.println(cdsStartIdxLeft+" "+cdsStartIdxRight+" "+cdsEndIdxLeft+" "+cdsEndIdxRight+" "+exons.size());
			
			Transcript t1 = null;
			Transcript t2 = null;
			
			if(cdsStartIdx > 0 && cdsStartIdx < exons.size()-1 && cdsStart > 500) {
				
				Transcript[] test = splitIntoTwo(cdsStartIdx,false);
				
				if(test[0].getNumberOfIntrons() > 1 && test[1].getNumberOfIntrons() > 1) {
					Transcript t1First = test[0];
					Transcript t2First = test[1];

					t2First.strand = this.strand;

					t1First.addStrandAndCDS(minProteinLength,true);
					t2First.addStrandAndCDS(minProteinLength,true);

					if(t1First.cdsEnd-t1First.cdsStart > 500 && t2First.cdsEnd-t2First.cdsStart > 500) {
						t1 = t1First;
						t2 = t2First;
					}
				}
				
			}
			
			if(cdsEndIdx > 0 && cdsEndIdx < exons.size()-1 && currPos-cdsEnd > 500) {
				
				Transcript[] test = splitIntoTwo(cdsEndIdx,false);
				
				if(test[0].getNumberOfIntrons() > 1 && test[1].getNumberOfIntrons() > 1) {
					Transcript t1Second = test[0];
					Transcript t2Second = test[1];

					t1Second.strand = this.strand;

					t1Second.addStrandAndCDS(minProteinLength,true);
					t2Second.addStrandAndCDS(minProteinLength,true);

					if(t1Second.cdsEnd-t1Second.cdsStart > 500 && t2Second.cdsEnd-t2Second.cdsStart > 500) {
						if(t1==null || t1Second.cdsEnd-t1Second.cdsStart + t2Second.cdsEnd-t2Second.cdsStart > t1.cdsEnd-t1.cdsStart+t2.cdsEnd-t2.cdsStart) {
							t1 = t1Second;
							t2 = t2Second;
						}
					}
				}
			}
			
			if(t1 != null) {
				res.add(new Transcript(t1.exons,t1.introns));
				res.add(new Transcript(t2.exons,t2.introns));
			}else {
				this.cdsStart = -1;
				this.cdsEnd = -1;
				this.strand = oldstrand;
				res.add(this);
			}
			
		}
		
		
		
		
		private Transcript[] splitIntoTwo(int possibleSplit, boolean isIntron) {
			LinkedList<Integer> exons1 = new LinkedList<Integer>();
			LinkedList<Integer> exons2 = new LinkedList<Integer>();
			LinkedList<Integer> introns1 = new LinkedList<Integer>();
			LinkedList<Integer> introns2 = new LinkedList<Integer>();
			
			
			
			for(int i=0;i<=possibleSplit&i<exons.size();i++) {
				exons1.add(this.exons.get(i));
			}
			
			for(int i=possibleSplit+(isIntron?1:0);i<exons.size();i++) {
				exons2.add(this.exons.get(i));
			}

			
			for(int i=0;i<this.introns.size();i++) {
				if(SplicingGraph.this.introns[this.introns.get(i)].right.regionStart <= contRegs[exons1.getLast()].regionStart) {
					introns1.add(this.introns.get(i));
				}else if(SplicingGraph.this.introns[this.introns.get(i)].left.regionEnd >= contRegs[exons2.getFirst()].regionEnd) {
					introns2.add(this.introns.get(i));
				}/*else {
					System.out.println("removed: "+(regionStart+SplicingGraph.this.introns[this.introns.get(i)].left.regionEnd)+" <-> "+(regionStart+SplicingGraph.this.introns[this.introns.get(i)].right.regionStart));
				}*/
			}
			
			Transcript t1 = new Transcript(exons1, introns1);
			t1.extraLeft = this.extraLeft;
			Transcript t2 = new Transcript(exons2, introns2);
			t2.extraRight = this.extraRight;
			
			return new Transcript[] {t1,t2};
		}
		
				
		
		public boolean testSplitBySpliceSitesAndORF2(LinkedList<Transcript> resL, int minProteinLength) {
			double[][] intronVotes = getIntronVotes();

			
			IntList idxs = new IntList();
			
			
			int i=0;
			int last = 0;
			while(i<intronVotes[0].length) {
				
				boolean didit = false;
				while(i<intronVotes[0].length && (intronVotes[0][i]>1 || (intronVotes[0][i]==0 && intronVotes[1][i]==0) )) {
					i++;
					didit = true;
				}
				
				if(i>last && i<intronVotes[0].length) {
					int j=i-1;
					while(j>0 && (intronVotes[0][j]==0 && intronVotes[1][j]==0)) {
						idxs.add(j-1);
						j--;
					}
					idxs.add(i-1);
					last = i;
				}
				
				while(i<intronVotes[1].length && (intronVotes[1][i]>1 || (intronVotes[0][i]==0 && intronVotes[1][i]==0))) {
					i++;
					didit = true;
				}
				if(i>last && i<intronVotes[1].length) {
					int j=i-1;
					while(j>0 && (intronVotes[0][j]==0 && intronVotes[1][j]==0)) {
						idxs.add(j-1);
						j--;
					}
					idxs.add(i-1);
					last = i;
				}
				if(!didit) {
					i++;
					last = i;
				}
			}
			
			
			int minIdx = 0;
			while( minIdx < intronVotes[0].length && intronVotes[0][minIdx] == 0 && intronVotes[1][minIdx] == 0 ) {
				minIdx++;
			}
			int maxIdx = intronVotes[0].length-1;
			while(maxIdx >= 0 && intronVotes[0][maxIdx] == 0 && intronVotes[1][maxIdx] == 0) {
				maxIdx--;
			}
			
			IntList idxs2 = new IntList();
			for(i=0;i<idxs.length();i++) {
				if(idxs.get(i)>= minIdx && idxs.get(i)<=maxIdx) {
					idxs2.add(idxs.get(i));
				}
			}
			
			
			
			if(idxs2.length() == 0) {
				resL.add(this);
				return false;
			}
			
			IntList pos = new IntList();
			
			
			
			for(i=0;i<idxs2.length();i++) {
				pos.add( regionStart + (contRegs[exons.get(idxs2.get(i))].regionStart+contRegs[exons.get(idxs2.get(i))].regionEnd)/2 );
			}
						
			pos.sort();
			
			
			int before = resL.size();
			
			testSplitByORF(resL, pos.toArray(), false,minProteinLength);
			

			return resL.size()>before;
			
		}
		
		

		private int getFirstIntronPredecessor() {
			int s1 = contRegs[exons.getFirst()].regionEnd;
			int x1 = exons.getFirst();
			int i=1;
			while(i<exons.size() && contRegs[exons.get(i)].regionStart==s1+1) {
				x1 = exons.get(i);
				s1 = contRegs[x1].regionEnd;
				i++;
			}
			if(i==exons.size()) {
				x1 = -1;
			}
			return x1;
		}
		
		private int getLastIntronSuccessor() {
			int e1 = contRegs[exons.getLast()].regionStart;
			int x1 = exons.getLast();
			int i=exons.size()-2;
			while(i>=0 && contRegs[exons.get(i)].regionEnd==e1-1) {
				x1 = exons.get(i);
				e1 = contRegs[x1].regionStart;
				i--;
			}
			if(i<0) {
				x1 = -1;
			}
			return x1;
		}

		public boolean equalsIntronChain(Transcript last) {
			if(last == null) {
				return false;
			}

			int x1 = getFirstIntronPredecessor();
			int x2 = last.getFirstIntronPredecessor();
			
			if(x1 == -1 && x2 == -1) {
				return true;
			}else if(x1 != x2){
				return false;
			}else{
				int y1 = getLastIntronSuccessor();
				int y2 = last.getLastIntronSuccessor();

				if(y1 != y2) {
					return false;
				}else {

					int i=0;
					while(i<exons.size() && exons.get(i) < x1) {
						i++;
					}
					int j=0;
					while(j<last.exons.size() && last.exons.get(j) < x2) {
						j++;
					}



					int k=0;
					while( i+k < exons.size() && j+k < last.exons.size() && exons.get(i+k) < y1 && last.exons.get(j+k) < y2  ) {
						if(exons.get(i+k) != last.exons.get(j+k)) {
							return false;
						}
						k++;
					}
					
					if(exons.get(i+k) != y1 || last.exons.get(j+k) != y2) {
						return false;
					}
					
					return true;

				}
			}
			
		}

		public boolean overlaps(Transcript last) {
			return this.getStart() >= last.getStart() && this.getStart() <= last.getEnd() ||
					this.getEnd() >= last.getStart() && this.getEnd() <= last.getEnd();
		}

		public double getMaximumNumberOfReads() {
			int num = 0;
			for(int i=0;i<exons.size();i++) {
				int idx = exons.get(i);
				if(contRegs[idx].getNumberOfReads()>num) {
					num = contRegs[idx].getNumberOfReads();
				}
			}
			return num;
		}

		public double getMaximumNumberOfCoveringIntronReads(LinkedList<Edge> intronEdges, int edgeRegionStart) {
			int start = getStart();
			int end = getEnd();
			
			int num = intronEdges.stream().filter( e -> e.getRelStart()+edgeRegionStart  < start && e.getRelEnd()+edgeRegionStart > end )
								.map(e -> e.getNumberOfReads())
								.reduce( 0, (a,b) -> Math.max(a, b));
			
			return num;
		}

		public int getCDSLength() {
			if(cdsStart<0) {
				return 0;
			}else {
				return cdsEnd-cdsStart+1;
			}
		}	
		
		public double getUTRFraction() {
			return 1.0 - getCDSLength()/(double)getLength();
		}
		
	}
	
	
	private class ContiguousRegion{
		
		private int regionStart;
		private int regionEnd;
		private LinkedList<Node> readNodes;
		private int[] readIdxs;
		private int numReads;
		
		private int[] firstReadIdxs;
		private int[] lastReadIdxs;
		
		private LinkedList<Intron> left;
		private LinkedList<Intron> right;
		
		public ContiguousRegion(int regionStart, int regionEnd, LinkedList<Node> readNodes) {
			super();
			this.regionStart = regionStart;
			this.regionEnd = regionEnd;
			this.readNodes = readNodes;
			
			this.numReads = 0;
			for(Node node : readNodes) {
				if(node.getNumberOfReads() > this.numReads) {
					this.numReads = node.getNumberOfReads();
				}
			}
			
			
			this.left = new LinkedList<SplicingGraph.Intron>();
			this.right = new LinkedList<SplicingGraph.Intron>();
		}
		
		public void setSortedReads(int[] readIdxs, int[] firstReadIdxs, int[] lastReadIdxs) {
			this.readIdxs = readIdxs;
			this.numReads = readIdxs.length;
			this.firstReadIdxs = firstReadIdxs;
			this.lastReadIdxs = lastReadIdxs;
		}
		
		public int getNumberOfReads() {
			return numReads;
		}
		
		public String toString() {
			return toString(SplicingGraph.this.regionStart);
		}
		
		public String toString(int offset) {
			return chromosome+"\tcontReg\texon\t"+(offset+regionStart)+"\t"+(offset+regionEnd)+"\t.\t+\t.\tID="+regionStart+"-"+regionEnd+"; left="+left.size()+"; right="+right.size()+"; cov="+getCoverage();
		}

		public double getCoverage() {
			double cov = 0;
			double n = 0;
			for(Node node : readNodes) {
				cov += node.getNumberOfReads();
				n++;
			}
			return cov/n;
		}
		
	}
	
	private class Intron implements Comparable {
		
		private ContiguousRegion left;
		private ContiguousRegion right;
		private Edge edge;
		private int[] readIdxs;
		
		public Intron(ContiguousRegion left, ContiguousRegion right, Edge edge) {
			this.left = left;
			this.right = right;
			this.edge = edge;
			
		}
		
		public double[] getFractions() {
			int[] leftReads = left.lastReadIdxs;
			int[] rightReads = right.firstReadIdxs;
			
			double leftIntersect = intersectSortedArrays(leftReads,readIdxs).length;
			double rightIntersect = intersectSortedArrays(rightReads,readIdxs).length;
			
			double nLeft = leftReads.length;
			double nRight = rightReads.length;
			
			return new double[] {leftIntersect/nLeft,rightIntersect/nRight};
			
			
		}
		

		public double getCoverage() {
			return edge.getNumberOfReads();
		}
		
		public String toString() {
			return toString(regionStart);
		}
		
		public String toString(int offset) {
			return chromosome+"\tintron\tintron\t"+(offset+left.regionEnd+1)+"\t"+(offset+right.regionStart-1)+"\t.\t+\t.\tID="+left.regionEnd+"-"+right.regionStart+"; left="+left.regionStart+"-"+left.regionEnd+"; right="+right.regionStart+"-"+right.regionEnd+"; cov="+getCoverage()+"; frac="+Arrays.toString( getFractions() );
		}
		
		public int getStart() {
			return left.regionEnd+1;
		}
		
		public int getEnd() {
			return right.regionStart-1;
		}
		
		public int getLength() {
			return right.regionStart-left.regionEnd-1;
		}

		@Override
		public int compareTo(Object o) {
			if(o instanceof Intron) {
				Intron i = (Intron)o;
				return Integer.compare(this.readIdxs.length, i.readIdxs.length);
			}else {
				return 1;
			}
		}
		
		
		public boolean isCanonical(int shortestIntronLength) {
			if(edge.getRelEnd()-edge.getRelStart()+1<shortestIntronLength) {
				return false;
			}else {
				char[] chr = Genome.genome.getChromosome(chromosome);
				int lastEnd = left.regionEnd;
				int start = right.regionStart;
				
				
				//String left = chr.substring(regionStart + lastEnd, regionStart + lastEnd+2);
				char left1 = chr[regionStart+lastEnd];
				char left2 = chr[regionStart+lastEnd+1];
				//String right = chr.substring(regionStart + start - 3, regionStart + start-1);
				char right1 = chr[regionStart+start-3];
				char right2 = chr[regionStart+start-2];
				
				boolean res = true;
				
				if(/*"GT".equals(left)*/
						'G' == left1 && 'T' == left2 || 
						'G' == left1 && 'C' == left2 || 
						'C' == left1 && 'T' == left2 || 
						'A' == left1 && 'T' == left2 ) {
					res &= true;
				}else {
					res &= false;
				}
				
				if(/*"AG".equals(right)*/
						'A' == right1 && 'G' == right2 || 
						'A' == right1 && 'C' == right2 || 
						'G' == right1 && 'C' == right2 ||
						'A' == right1 && 'C' == right2 ||
						'A' == right1 && 'T' == right2) {
					res &= true;
				}else {
					res &= false;
				}
				
				return res;
				
			}
		}
		
		public char getStrandBySpliceSites(int shortestIntronLength) {
			if(edge.getRelEnd()-edge.getRelStart()+1<shortestIntronLength) {
				return '.';
			}else {
				char[] chr = Genome.genome.getChromosome(chromosome);
				int lastEnd = left.regionEnd;
				int start = right.regionStart;
				
				//String left = chr.substring(regionStart + lastEnd, regionStart + lastEnd+2);
				char left1 = chr[regionStart+lastEnd];
				char left2 = chr[regionStart+lastEnd+1];
				//String right = chr.substring(regionStart + start - 3, regionStart + start-1);
				char right1 = chr[regionStart+start-3];
				char right2 = chr[regionStart+start-2];

				if(/*"GT".equals(left) && "AG".equals(right)*/
						'G' == left1 && 'T' == left2 && 'A' == right1 && 'G' == right2) {
					return '+';
				}else if(/*"CT".equals(left) && "AC".equals(right)*/
						'C' == left1 && 'T' == left2 && 'A' == right1 && 'C' == right2) {
					return '-';
				}else {
					return '.';
				}
				
			}
		}
		
	}
	
	
	
	public class TranscriptLasso extends DifferentiableFunction{

		private int[][] assign;
		private double[] totalLen;
		private double[] target;
		private double[] vals;
		private double lambda1, lambda2;
		
		public TranscriptLasso(int[][] assign, double[] target, double lambda1, double lambda2) {
			super();
			this.assign = assign;
			this.totalLen = new double[assign.length];
			for(int i=0;i<assign.length;i++) {
				for(int j=0;j<assign[i].length;j++) {
					totalLen[i] += contRegs[ assign[i][j] ].regionEnd - contRegs[ assign[i][j] ].regionStart + 1; 
				}
			}
			this.target = target;
			this.vals = new double[target.length];
			this.lambda1 = lambda1;
			this.lambda2 = lambda2;
		}

		@Override
		public double evaluateFunction(double[] weights) throws DimensionException, EvaluationException {
			Arrays.fill(this.vals, 0.0);
			
			for(int i=0;i<assign.length;i++) {
				for(int j=0;j<assign[i].length;j++) {
					vals[ assign[i][j] ] += Math.exp( weights[i] ) / totalLen[i] * ( contRegs[assign[i][j]].regionEnd-contRegs[assign[i][j]].regionStart+1);
				}
			}
			
			double val = 0.0;
			for(int i=0;i<target.length;i++) {
				double temp = vals[i]*Math.exp( weights[ weights.length-1 ] ) - target[i];
				val += temp*temp;
			}
			
			for(int i=0;i<weights.length;i++) {
				val += lambda1*Math.abs(Math.exp( weights[i] ));
				
				val += lambda2*Math.exp( 2.0*weights[i] );
			}
			
			return val;
			
		}

		@Override
		public int getDimensionOfScope() {
			return assign.length+1;
		}

		@Override
		public double[] evaluateGradientOfFunction(double[] weights) throws DimensionException, EvaluationException {
			
			double[] grad = new double[weights.length];			
			
			Arrays.fill(this.vals, 0.0);
			
			for(int i=0;i<assign.length;i++) {
				for(int j=0;j<assign[i].length;j++) {
					vals[ assign[i][j] ] += Math.exp( weights[i] )/ totalLen[i] * ( contRegs[assign[i][j]].regionEnd-contRegs[assign[i][j]].regionStart+1);
				}
			}
			
			System.out.println(weights[weights.length-1]+" "+Arrays.toString(vals));
			
			double[] temps = new double[target.length];
			for(int i=0;i<target.length;i++) {
				temps[i] = 2.0 * ( vals[i]*Math.exp( weights[ weights.length-1 ] ) - target[i] );
				
			}
			
			for(int i=0;i<assign.length;i++) {
				for(int j=0;j<assign[i].length;j++) {
					grad[i] += temps[assign[i][j]] * Math.exp(weights[i]);
					grad[grad.length-1] += temps[assign[i][j]] * vals[assign[i][j]] * Math.exp( weights[ weights.length-1 ] );
				}
			}
			
			for(int i=0;i<weights.length;i++) {
				grad[i] += Math.exp(weights[i])*lambda1;
				grad[i] += 2*Math.exp(2.0*weights[i])*lambda2;
			}
			
			return grad;
			
			
		}
		
		
		
		
		
	}

	
	public LinkedList<Transcript> testSplitByORF(LinkedList<Transcript> list, int[] proposedSplits, int minProteinLength){
		LinkedList<Transcript> res = new LinkedList<SplicingGraph.Transcript>();
		boolean split = false;
		for(Transcript t : list) {
			char tempStrand = t.getStrand();
			t.addStrandAndCDS(minProteinLength,true);

			double utrFrac = 1.0 - (t.cdsEnd-t.cdsStart)/(double)t.getLength();
			t.setStrand(tempStrand);
			t.cdsStart = t.cdsEnd = -1;
						
			if(t.getNumberOfIntrons() < 1 && utrFrac < 0.5) {
				res.add(t);
			}else {
				int s = res.size();
				t.testSplitByORF(res,proposedSplits, true,minProteinLength);
				split |= res.size()-s>1;
			}
		}
	
		return res;
		
	}
	

	public LinkedList<Transcript> testSplitByORF(LinkedList<Transcript> list, int minProteinLength){
		LinkedList<Transcript> res = new LinkedList<SplicingGraph.Transcript>();
		boolean split = false;
		for(Transcript t : list) {
			if(t.getNumberOfIntrons() < 1) {
				res.add(t);
			}else {
				int s = res.size();
				t.testSplitByORF2(res,minProteinLength);
				split |= res.size()-s>1;
			}
		}
			
		return res;
		
	}
	
	
	public LinkedList<Transcript> testSplitBySpliceSitesAndORF(LinkedList<Transcript> list, int minProteinLength) {
		LinkedList<Transcript> res = new LinkedList<SplicingGraph.Transcript>();
		boolean split = false;
		for(Transcript t : list) {
			if(t.getNumberOfIntrons() == 0) {
				res.add(t);
			}else {
				//int s = res.size();
				boolean temp = t.testSplitBySpliceSitesAndORF2(res,minProteinLength);
				split |= temp;
				//split |= res.size()-s>1;
			}
		}
		
		return res;
		
	}
	
	private void makeUnique(LinkedList<Transcript> res) {
		Collections.sort(res,new Comparator<Transcript>() {

			@Override
			public int compare(Transcript o1, Transcript o2) {
				int s1 = o1.getFirstIntronPredecessor();
				int s2 = o2.getFirstIntronPredecessor();
				int cmp = Integer.compare(s1, s2);
				if(cmp == 0) {
					int e1 = o1.getLastIntronSuccessor();
					int e2 = o2.getLastIntronSuccessor();
					cmp = Integer.compare(e1, e2);
					if(cmp == 0) {
						
						
						int x1 = o1.getFirstIntronPredecessor();
						int x2 = o2.getFirstIntronPredecessor();
						
						int i=0;
						while(i<o1.exons.size() && o1.exons.get(i) < x1) {
							i++;
						}
						int j=0;
						while(j<o2.exons.size() && o2.exons.get(j) < x2) {
							j++;
						}
						
						int k=0;
						while( i+k < o1.exons.size() && j+k < o2.exons.size() ) {
							cmp = Integer.compare(o1.exons.get(i+k), o2.exons.get(j+k));
							if(cmp != 0) {
								return cmp;
							}
							k++;
						}
						
					}
				}
				return cmp;
			}
			
		});
		
		
		IntList toRemove = new IntList();
		int i=0;
		Transcript last = null;
		int lastIdx = -1;
		for(Transcript t : res) {
			
			if(t.equalsIntronChain(last) && t.overlaps(last)) {
				
				if(last.getEnd()-last.getStart() > t.getEnd() - t.getStart()) {
					toRemove.add(i);
				}else {
					if(toRemove.length()== 0 || toRemove.get(toRemove.length()-1) != lastIdx) {
						toRemove.add(lastIdx);
					}
					last = t;
					lastIdx = i;
				}
			}else {
				last = t;
				lastIdx = i;
			}
			i++;
		}
		toRemove.sort();

		for(i=toRemove.length()-1;i>=0;i--) {
			res.remove(toRemove.get(i));
		}
		
	}

	
	public LinkedList<Transcript>[] splitByLocationAndMakeUnique(LinkedList<Transcript> list, int minProteinLength, Stranded stranded) {
		Collections.sort(list);

		int endPlus = 0;
		int endMinus = 0;
		int endOther = 0;
		
		LinkedList<LinkedList<Transcript>> lplus = null;
		LinkedList<LinkedList<Transcript>> lminus = null;
		LinkedList<LinkedList<Transcript>> lother = null;
		
		
		for(Transcript t : list) {
			
			char strand = t.getStrand();
			if(strand == '.') {
				strand = t.getStrandByIntrons();
			}
			
			if(strand == '+') {
				
				if(lplus == null) {
					lplus = new LinkedList<LinkedList<Transcript>>();
				}
				if(t.getStart()>endPlus) {
					lplus.add( new LinkedList<SplicingGraph.Transcript>() );
				}
				lplus.getLast().add(t);
				if(t.getEnd() > endPlus) {
					endPlus = t.getEnd();
				}
				
			}else if(strand == '-') {
				
				if(lminus == null) {
					lminus = new LinkedList<LinkedList<Transcript>>();
				}
				if(t.getStart()>endMinus) {
					lminus.add( new LinkedList<SplicingGraph.Transcript>() );
				}
				lminus.getLast().add(t);
				if(t.getEnd() > endMinus) {
					endMinus = t.getEnd();
				}
			}else {
				int diffplus = t.getStart()-endPlus;
				int diffminus = t.getStart()-endMinus;
				
				if(stranded == Stranded.FR_UNSTRANDED && ( (lplus!= null && diffplus < 0) || (lminus != null && diffminus < 0) ) ) {
					if(diffplus < diffminus) {
						lplus.getLast().add(t);
					}else {
						lminus.getLast().add(t);
					}
				}else {
					if(lother == null) {
						lother = new LinkedList<LinkedList<Transcript>>();
					}
					if(t.getStart()>endOther) {
						lother.add( new LinkedList<SplicingGraph.Transcript>() );
					}
					lother.getLast().add(t);
					if(t.getEnd() > endOther) {
						endOther = t.getEnd();
					}
				}
			}
			
			t.addStrandAndCDS(minProteinLength,false);
		}
		
		
		if(lother != null && lother.size() == 1) {
			if(lplus == null && lminus != null && lminus.size() == 1) {
				lminus.getFirst().addAll(lother.getFirst());
				lother = null;
			}else if(lplus != null && lminus == null && lplus.size() == 1) {
				lplus.getFirst().addAll(lother.getFirst());
				lother = null;
			}
		}
		
		LinkedList<LinkedList<Transcript>> res = new LinkedList<LinkedList<Transcript>>(); 
		
		if(lother != null) {
			res.addAll(lother);
		}
		if(lplus != null) {
			res.addAll(lplus);
		}
		if(lminus != null) {
			res.addAll(lminus);
		}
		
		
		for(LinkedList<Transcript> li : res) {
			makeUnique(li);
		}
		
		
		
		return res.toArray(new LinkedList[0]);
		
	}
	
	public LinkedList<Transcript>[] splitByStrandAndMakeUnique(LinkedList<Transcript> sList, int minProteinLength) {
		LinkedList<Transcript>[] sList2 = new LinkedList[] {
				new LinkedList<Transcript>(),
				new LinkedList<Transcript>(),
				new LinkedList<Transcript>()
		};
		for(Transcript t : sList) {
			
			t.addStrandAndCDS(minProteinLength,false);
			if(t.getStrand() == '+') {
				sList2[1].add(t);
			}else if(t.getStrand() == '-') {
				sList2[2].add(t);
			}else {
				sList2[0].add(t);
			}
			
		}
		if(sList2[1].size() == 0 && sList2[2].size() > 0) {
			sList2[2].addAll(sList2[0]);
			sList2[0].clear();
		}else if(sList2[1].size() > 0 && sList2[2].size() == 0) {
			sList2[1].addAll(sList2[0]);
			sList2[0].clear();
		}
		
		makeUnique(sList2[0]);
		makeUnique(sList2[1]);
		makeUnique(sList2[2]);
		
		return sList2;
		
	}
	
	
	
	
}
