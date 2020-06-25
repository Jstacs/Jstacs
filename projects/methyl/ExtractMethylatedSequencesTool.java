/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package projects.methyl;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.zip.GZIPInputStream;

import de.jstacs.DataType;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.DataColumnParameter;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;

public class ExtractMethylatedSequencesTool implements JstacsTool {

	
	public static void main(String[] args) throws Exception{
		CLI cli = new CLI(new ExtractMethylatedSequencesTool());
		
		cli.run(args);
	}
	
	private static class Peak{
		
		public int center;
		public double stat;
		
		
		public Peak(int center, double stat) {
			super();
			this.center = center;
			this.stat = stat;
		}


		public int getCenter() {
			return center;
		}


		public double getStat() {
			return stat;
		}
		
		
		
	}
	
	
	public ExtractMethylatedSequencesTool() {
		// TODO Auto-generated constructor stub
	}

	@Override
	public ToolParameterSet getToolParameters() {
		
		LinkedList<Parameter> parameters = new LinkedList<Parameter>();
		
		parameters.add(new FileParameter("Genome", "The FastA containing all chromosome sequences, may be gzipped", "fa,fas,fasta,fa.gz,fas.gz,fasta.gz", true));
		
		FileParameter peaks = new FileParameter("Peaks", "The file containing the peaks in tabular format", "bed,gff,gff3,narrowPeak,gtf,tabular", true);
		
		parameters.add(peaks);
		
		try{
		
			parameters.add(new DataColumnParameter(peaks.getName(), "Chromosome column", "The column of the peaks file containing the chromosome", true,1));

			parameters.add(new DataColumnParameter(peaks.getName(), "Start column", "The column of the peaks file containing the start position relative to the chromsome start", true,2));

			SelectionParameter sp2 = new SelectionParameter(DataType.PARAMETERSET,new String[]{
					"Peak center",
					"End of peak"
			},new Object[]{
					new SimpleParameterSet(new DataColumnParameter(peaks.getName(), "Center column", "The column of the peaks file containing the peak center relative to the start position", true)),
					new SimpleParameterSet(new DataColumnParameter(peaks.getName(), "End column", "The column of the peaks file containing the end position relative to the chromsome start", true,3))
			},"Peak position", "The kind how the peak is specified", true );
			sp2.setDefault("End of peak");
			
			parameters.add(sp2);

			
			parameters.add(new SimpleParameter(DataType.INT, "Width", "The fixed width of all extracted regions",true, new NumberValidator<Integer>(1, 10000),1000));
			
			parameters.add(new DataColumnParameter(peaks.getName(),"Statistics column","The column of the peaks file containing the peak statistic or a similar measure of confidence",true,7));
			
			
		}catch(Exception doesnothappen){throw new RuntimeException();}
		
		return new ToolParameterSet(getShortName(),parameters);
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		
		progress.setLast(1.0);
		progress.setCurrent(0.0);
		
		String genome = ((FileParameter)parameters.getParameterAt(0)).getFileContents().getFilename();
		
		String peaks = ((FileParameter)parameters.getParameterAt(1)).getFileContents().getFilename();
		
		int chromcol = (Integer) parameters.getParameterAt(2).getValue();
		
		int startcol = (Integer) parameters.getParameterAt(3).getValue();
		
		SelectionParameter sp = (SelectionParameter) parameters.getParameterAt(4);
		
		boolean isCenter=false;
		int seccol = 0;
		
		if(sp.getSelected() == 0){
			isCenter = true;
			seccol = (Integer) ((ParameterSet)sp.getValue()).getParameterAt(0).getValue();
		}else{
			isCenter = false;
			seccol = (Integer) ((ParameterSet)sp.getValue()).getParameterAt(0).getValue();
		}
		
		int width = (Integer) parameters.getParameterAt(5).getValue();
		width /=2;
		
		int statcol = (Integer) parameters.getParameterAt(6).getValue();
		
		
		HashMap<String, LinkedList<Peak>> map = new HashMap<String, LinkedList<Peak>>();
		
		BufferedReader read = new BufferedReader(new FileReader(peaks));
		
		String str = null;
		while( (str = read.readLine()) != null ){
			String[] parts = str.split("\t");
			String chrom = parts[chromcol-1];
			int start = Integer.parseInt(parts[startcol-1]);
			int sec = Integer.parseInt(parts[seccol-1]);
			double stat = Double.parseDouble(parts[statcol-1]);
			
			int center;
			if(isCenter){
				center = start+sec;
			}else{
				center = (start+sec)/2;
			}
			
			if(!map.containsKey(chrom)){
				map.put(chrom, new LinkedList<ExtractMethylatedSequencesTool.Peak>());
			}
			
			map.get(chrom).add(new Peak(center,stat));
			
		}
		
		progress.setCurrent(0.1);
		
		read.close();
		
		if(genome.toLowerCase().endsWith(".gz")) {
			read = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(genome))));
		}else {
			read = new BufferedReader(new FileReader(genome));
		}
		String chrom = null;
		
		
		StringBuffer currChrom = new StringBuffer();
		
		StringBuffer res = new StringBuffer();
		int i=1;
		int num = map.keySet().size();
		int nSkipped = 0;
		while( true ){
			str = read.readLine();
			if(str ==null || str.startsWith(">")){
				if(chrom != null){
					//System.out.println(chrom);
					LinkedList<Peak> myPeaks = map.remove(chrom);
					if(myPeaks != null){
						
						Iterator<Peak> it = myPeaks.iterator();
						while(it.hasNext()){
							Peak peak = it.next();
							int start = peak.getCenter()-width-1;
							int end = peak.getCenter()+width-1;
							if(start>=0 && end <= currChrom.length()){
								String sub = currChrom.substring(start, end);
								if(sub.matches("^[ACGTHMacgthm]+$")){
									res.append(">chrom: "+chrom+"; center: "+peak.getCenter()+"; peak: "+(width+1)+"; signal: "+peak.getStat()+"\n");
									res.append(sub+"\n");
								}else{
									nSkipped++;
									//protocol.appendWarning("Peak at "+chrom+":"+start+"-"+end+" skipped because of ambiguous nucleotides.\n");
									//System.out.println(">chrom: "+chrom+"; center: "+peak.getCenter()+"; peak: "+(width+1)+"; signal: "+peak.getStat()+"\n");
									//System.out.println(sub+"\n");
									
								}
							}else{
								protocol.appendWarning("Peak at "+chrom+":"+start+"-"+end+" spans outside chromsome "+chrom+" of length "+currChrom.length()+".\n");
							}
						}
						
						progress.setCurrent(0.1 + 0.8*i/(double)(num+1));

						i++;
					}
					if(nSkipped > 0) {
						protocol.appendWarning(nSkipped+" peak"+(nSkipped > 1 ? "s" : "")+" on chromosome "+chrom+" skipped because of ambiguous nucleotides.\n");
					}
				}
				
				currChrom.delete(0, currChrom.length());
				nSkipped = 0;
				if(str == null){
					break;
				}else{
					int space = str.indexOf(" ");
					chrom = str.substring(1,space < 0 ? str.length() : space).trim();
					if(!map.containsKey(chrom)){
						protocol.appendWarning("No peaks on chromosome "+chrom+".\n");
					}
				}
			}else{
				currChrom.append(str.trim());
			}
		}
		
		Iterator<String> missed = map.keySet().iterator();
		
		while( missed.hasNext() ){
			protocol.appendWarning("No sequence for "+missed.next()+".\n");
		}
		
		TextResult extracted = new TextResult("Extracted sequences", "The sequences under the peaks", new FileParameter.FileRepresentation("", res.toString()), "fasta", "ExtractSequences", null, true);
		
		protocol.append("Extraction finished\n");
		
		progress.setCurrent(1.0);
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(extracted), parameters, getToolName(), new Date(System.currentTimeMillis()));
	}

	@Override
	public String getToolName() {
		return "Data Extractor";
	}

	@Override
	public String getToolVersion() {
		return "1.0";
	}

	@Override
	public String getShortName() {
		return "extract";
	}

	@Override
	public String getDescription() {
		return "extracts data from a genome and peak coordinates as required for Dimont";
	}

	@Override
	public String getHelpText() {
		return "**Data Extractor** prepares an annotated FastA file as required by Dimont from a genome (in FastA format, including methylated variants) and a tabular file (e.g., BED, GTF, narrowPeak,...). "
				+ "The regions specified in the tabular file are used to determine the center of the extracted sequences. All extracted sequences have the same length as specified by parameter \"Width\".\n" + 
				"\n" + 
				"In case of ChIP data, the center position could for instance be the peak summit.\n" + 
				"An annotated FastA file for ChIP-seq data comprising sequences of length 100 centered around the peak summit might look like::\n" + 
				"	\n" + 
				"	> peak: 50; signal: 515\n" + 
				"	ggccatgtgtatttttttaaatttccac...\n" + 
				"	> peak: 50; signal: 199\n" + 
				"	GGTCCCCTGGGAGGATGGGGACGTGCTG...\n" + 
				"	...\n" + 
				"\n" + 
				"where the center is given as 50 for the first two sequences, and the confidence amounts to 515 and 199, respectively.\n" +
				"\n" + 
				"\n" + 
				"If you experience problems using Data Extractor, please contact_ us.\n" + 
				"\n" + 
				".. _contact: mailto:grau@informatik.uni-halle.de";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[]{new ResultEntry(TextResult.class, "fasta", "Extracted sequences")};
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
