package projects.dimont;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.DataSetResult;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.DataColumnParameter;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolResult;
import projects.motifComp.FindPWMsAndClusters;

public class ExtractSequencesTool implements JstacsTool {

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
	
	
	public ExtractSequencesTool() {
		// TODO Auto-generated constructor stub
	}

	@Override
	public ParameterSet getToolParameters() {
		
		LinkedList<Parameter> parameters = new LinkedList<Parameter>();
		
		parameters.add(new FileParameter("Genome", "The FastA containing all chromosome sequences", "fa,fas,fasta", true));
		
		FileParameter peaks = new FileParameter("Peaks", "The file containing the peaks in tabular format", "bed,gff,gff3,narrowPeak,gtf,tabular", true);
		
		parameters.add(peaks);
		
		try{
		
			parameters.add(new DataColumnParameter(peaks.getName(), "Chromosome column", "The column of the peaks file containing the chromosome", true));

			parameters.add(new DataColumnParameter(peaks.getName(), "Start column", "The column of the peaks file containing the start position relative to the chromsome start", true));

			parameters.add(new SelectionParameter(DataType.PARAMETERSET,new String[]{
					"Peak center",
					"End of peak"
			},new Object[]{
					new SimpleParameterSet(new DataColumnParameter(peaks.getName(), "Center column", "The column of the peaks file containing the peak center relative to the start position", true)),
					new SimpleParameterSet(new DataColumnParameter(peaks.getName(), "End column", "The column of the peaks file containing the end position relative to the chromsome start", true))
			},"Peak position", "The kind how the peak is specified", true ));

			parameters.add(new SimpleParameter(DataType.INT, "Width", "The fixed width of all extracted regions",true, new NumberValidator<Integer>(1, 10000),1000));
			
			parameters.add(new DataColumnParameter(peaks.getName(),"Statistics column","The column of the peaks file containing the peak statistic or a similar measure of confidence",true));
			
			
		}catch(Exception doesnothappen){throw new RuntimeException();}
		
		return new SimpleParameterSet(parameters.toArray(new Parameter[0]));
		
	}

	@Override
	public ToolResult run(ParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		
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
				map.put(chrom, new LinkedList<ExtractSequencesTool.Peak>());
			}
			
			map.get(chrom).add(new Peak(center,stat));
			
		}
		
		progress.setCurrent(0.1);
		
		read.close();
		
		read = new BufferedReader(new FileReader(genome));
		
		String chrom = null;
		
		
		StringBuffer currChrom = new StringBuffer();
		
		StringBuffer res = new StringBuffer();
		int i=1;
		int num = map.keySet().size();
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
								if(sub.matches("^[ACGTacgt]+$")){
									res.append(">chrom: "+chrom+"; center: "+peak.getCenter()+"; peak: "+(width+1)+"; signal: "+peak.getStat()+"\n");
									res.append(sub+"\n");
								}else{
									protocol.appendWarning("Peak at "+chrom+":"+start+"-"+end+" skipped because of ambiguous nucleotides.\n");
								}
							}else{
								protocol.appendWarning("Peak at "+chrom+":"+start+"-"+end+" spans outside chromsome "+chrom+" of length "+currChrom.length()+".\n");
							}
						}
						
						progress.setCurrent(0.1 + 0.8*i/(double)(num+1));

						i++;
					}
				}
				currChrom.delete(0, currChrom.length());
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
		try {
			return FileManager.readInputStream( FindPWMsAndClusters.class.getClassLoader().getResourceAsStream( "projects/dimont/helpExtractor.txt" ) ).toString();
		} catch ( Exception e ) {
			e.printStackTrace();
			return "";
		}
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[]{new ResultEntry(TextResult.class, "fasta", "Extracted sequences")};
	}

}
