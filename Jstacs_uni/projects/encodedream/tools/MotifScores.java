package projects.encodedream.tools;
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
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList;
import java.util.Date;
import java.util.LinkedList;
import java.util.zip.GZIPOutputStream;

import de.jstacs.DataType;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.AbstractSelectionParameter.InconsistentCollectionException;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.sequenceScores.QuickScanningSequenceScore;
import de.jstacs.sequenceScores.statisticalModels.trainable.PFMWrapperTrainSM;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.PFMComparator;
import projects.dimont.ThresholdedStrandChIPper;
import projects.encodedream.AggregateMotifProfiles;
import projects.encodedream.QuickMotifProfileTool;
import projects.encodedream.SlowMotifProfileTool;

public class MotifScores implements JstacsTool {

	@Override
	public ToolParameterSet getToolParameters() {
		
		LinkedList<Parameter> pars = new LinkedList<>();
		
		try {
			pars.add(new SelectionParameter(DataType.PARAMETERSET,new String[]{"Dimont","HOCOMOCO","Jaspar"},
					new ParameterSet[]{
							new SimpleParameterSet( new FileParameter("Dimont motif", "Dimont motif model description", "xml", true) ),
							new SimpleParameterSet( new FileParameter("HOCOMOCO PWM", "PWM from the HOCOMOCO database", "txt,pwm", true) ),
							new SimpleParameterSet( new FileParameter("Jaspar PFM", "PFM in Jaspar format", "txt", true) )
					},"Motif model","The motif model in Dimont, HOCOMOCO, or Jaspar format",true));
		} catch (InconsistentCollectionException | IllegalValueException | DatatypeNotValidException e1) {
			e1.printStackTrace();
		}
		
		//pars.add(new FileParameter("Motif model", "Dimont motif model description", "xml", true));
		pars.add(new FileParameter("Genome","Genome as FastA file","fa,fas,fasta",true));
		pars.add(new FileParameter("FAI of genome", "FastA index file of the genome", "fai", true));
		try {
			pars.add(new SimpleParameter(DataType.INT, "Bin width", "The width of the genomic bins considered", true));
			pars.add(new SimpleParameter(DataType.BOOLEAN, "Low-memory mode", "Use slower mode with a smaller memory footprint", true, false));
		} catch (DatatypeNotValidException | IllegalValueException e) {
			e.printStackTrace();
		}
		
		return new ToolParameterSet(getShortName(),pars.toArray(new Parameter[0]));
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		QuickScanningSequenceScore model = null;
		if(((SelectionParameter)parameters.getParameterAt(0)).getSelected()==0){
			String slimfile = (String) ((ParameterSet)parameters.getParameterAt(0).getValue()).getParameterAt(0).getValue();
			GenDisMixClassifier cl = new GenDisMixClassifier(FileManager.readFile(slimfile));
			ThresholdedStrandChIPper fg =(ThresholdedStrandChIPper) cl.getDifferentiableSequenceScore(0);
			model = (QuickScanningSequenceScore) fg.getFunction(0);
			
		}else if(((SelectionParameter)parameters.getParameterAt(0)).getSelected()==1){
			
			String hocofile = (String) ((ParameterSet)parameters.getParameterAt(0).getValue()).getParameterAt(0).getValue();
			
			BufferedReader reader = new BufferedReader(new FileReader(hocofile));
			
			reader.readLine();
			LinkedList<double[]> lines = new LinkedList<>();
			String str = null;
			while( (str = reader.readLine()) != null ){
				String[] parts = str.split("\t");
				double[] temp = new double[parts.length];
				for(int j=0;j<parts.length;j++){
					temp[j] = Double.parseDouble( parts[j] );
				}
				lines.add(temp);
			}
			reader.close();
			
			double[][] pssm = lines.toArray(new double[0][]);
			
			model = new PFMWrapperTrainSM(DNAAlphabetContainer.SINGLETON, "", pssm);
			
			
		}else{
			
			String jasfile = (String) ((ParameterSet)parameters.getParameterAt(0).getValue()).getParameterAt(0).getValue();
			
			ArrayList<SimpleEntry<String, double[][]>> pwms = PFMComparator.readPFMsFromJasparFastA( new BufferedReader( new FileReader(jasfile) ) );
			
			model = new PFMWrapperTrainSM(DNAAlphabetContainer.SINGLETON,null,pwms.get(0).getValue(),0.0);
			
		}
		
		
		final QuickScanningSequenceScore model2 = model;
		
		String genome = ((FileParameter)parameters.getParameterAt(1)).getFileContents().getFilename();
		String faiFile = ((FileParameter)parameters.getParameterAt(2)).getFileContents().getFilename();
		int bin = (int) parameters.getParameterAt(3).getValue();
		boolean lowmem = (boolean) parameters.getParameterAt(4).getValue();
		
		
		
		PipedInputStream in = new PipedInputStream();
		PipedOutputStream out = new PipedOutputStream(in);
		
		new Thread(()->{
			try {
				if(lowmem){
					SlowMotifProfileTool mot = new SlowMotifProfileTool();
					mot.run(model2, genome, threads, new BufferedOutputStream(out));
					out.close();
				}else{
					QuickMotifProfileTool mot = new QuickMotifProfileTool();
					mot.run(model2, genome, threads, new BufferedOutputStream(out));
					out.close();
				}
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		}).start();
		
		File outfile = File.createTempFile("motif", ".temp.gz", new File("."));
		outfile.deleteOnExit();
		
		PrintStream out2 = new PrintStream(new GZIPOutputStream(new FileOutputStream(outfile)));
		
		AggregateMotifProfiles.run(new BufferedReader(new InputStreamReader(in)), out2, bin,faiFile);
		
		in.close();
		out2.close();
		
		TextResult tr = new TextResult("Motif scores", "Features computed from the profile of motif scores", new FileParameter.FileRepresentation(outfile.getAbsolutePath()), "tsv.gz", getToolName(), null, true);
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );
		
		
	}

	@Override
	public String getToolName() {
		return "Motif scores";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "motif";
	}

	@Override
	public String getDescription() {
		return "computes motif-based features";
	}

	@Override
	public String getHelpText() {
		return "**Motif scores** computes features based on motif scores of a given motif model scanning sub-sequences along the genome. Motif scores are aggregated "
				+ "in bins of the specified width as maximum score and log of the average exponential score (i.e., average log-likelihood in case of statistical models). The motif model may be provided "
				+ "as PWMs in HOCOMOCO or PFMs in Jaspar format, or as Dimont motif models in XML format. For more complex motif models like Slim models, the current implementation uses several indexes to "
				+ "speed-up the scanning process. However, computation of these indexes is rather memory-consuming and often not reasonable for simple PWM models. Hence, a low-memory variant of the tool is available, "
				+ "which is typically only slightly slower for PWM models but substantially slower for Slim models. Output is provided as a gzipped file *Motif_scores.tsv.gz* containing columns chromosome, start position, "
				+ "maximum and average score. This output file together with a protocol of the tool run is saved to the specified output directory.";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}

}
