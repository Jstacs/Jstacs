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
import java.io.File;
import java.util.Date;
import java.util.HashSet;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.io.FileManager;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import projects.encodedream.ClassifiersWithInfo;
import projects.encodedream.FeatureReader;
import projects.encodedream.IterativeTraining;

public class IterativeTrainingTool implements JstacsTool {

	@Override
	public ToolParameterSet getToolParameters() {
		LinkedList<Parameter> pars = new LinkedList<>();
		
		pars.add( new FileParameter("Accessibility", "File containing accessibility features", "tsv.gz", true));
		try {
			pars.add( new ParameterSetContainer( new ExpandableParameterSet(new SimpleParameterSet(
						new FileParameter("Motif", "File containing motif features", "tsv.gz", true)
					), "Motif features", "File(s) containing the motif features")) );
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		
		pars.add( new FileParameter("Labels", "File containing the labels", "tsv.gz", true));
		
		pars.add(new FileParameter("FAI of genome", "FastA index file of the genome", "fai", true));
		
		try {
			pars.add(new SimpleParameter(DataType.INT, "Bin width", "The width of the genomic bins", true, new NumberValidator<Integer>(1, 1000), 50));
			pars.add(new SimpleParameter(DataType.INT, "Number of bins", "The number of adjacent bins", true, new NumberValidator<Integer>(1, 20), 5));
			pars.add(new SimpleParameter(DataType.INT, "Aggregation: bins before", "The number of bins before the current one considered in the aggregation", true, new NumberValidator<Integer>(1, 20), 1));
			pars.add(new SimpleParameter(DataType.INT, "Aggregation: bins after", "The number of bins after the current one considered in the aggregation", true, new NumberValidator<Integer>(1, 20), 4));
			pars.add(new SimpleParameter(DataType.INT, "Iterations", "The number of iterations of the interative training", true, new NumberValidator<Integer>(1, 20), 5));
			
			pars.add(new SimpleParameter(DataType.STRING, "Training chromosomes", "Training chromosomes, separated by commas", false));
			pars.add(new SimpleParameter(DataType.STRING, "Iterative training chromosomes", "Chromosomes with predictions in iterative training, separated by commas", false));
			
			pars.add(new SimpleParameter(DataType.DOUBLE, "Percentile", "Percentile of the prediction scores of positives used as threshold in iterative training", true, new NumberValidator<Double>(0.0, 1.0), 0.01));
		} catch (ParameterException e) {
			e.printStackTrace();
		}
		
		
		return new ToolParameterSet(getShortName(),pars.toArray(new Parameter[0]));
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		int numBins = (int) parameters.getParameterAt(5).getValue();
		int bin = (int) parameters.getParameterAt(4).getValue();
		String labelsFile = (String) parameters.getParameterAt(2).getValue();
		String dnaseFile = (String) parameters.getParameterAt(0).getValue();
		
		int num = ((ExpandableParameterSet) parameters.getParameterAt(1).getValue() ).getNumberOfParameters();
		String[] motifFiles = new String[num];
		for(int i=0;i<num;i++){
			motifFiles[i] = (String) ((ParameterSet)((ExpandableParameterSet) parameters.getParameterAt(1).getValue() ).getParameterAt(i).getValue()).getParameterAt(0).getValue();
		}
		
		String faiFile = (String) parameters.getParameterAt(3).getValue();
		int numBefore = (int) parameters.getParameterAt(6).getValue();
		int numAfter = (int) parameters.getParameterAt(7).getValue();
		int iterations = (int) parameters.getParameterAt(8).getValue();
		HashSet<String> trainChroms = null;
		if(((String)parameters.getParameterAt(9).getValue()) != null && ((String)parameters.getParameterAt(9).getValue()).length()>0){
			trainChroms = new HashSet<>();
			String[] parts = ((String)parameters.getParameterAt(9).getValue()).split(",");
			for(int i=0;i<parts.length;i++){
				trainChroms.add(parts[i]);
			}
		}
		LinkedList<String> itChroms = new LinkedList<>(trainChroms);
		if(((String)parameters.getParameterAt(10).getValue()) != null && ((String)parameters.getParameterAt(10).getValue()).length()>0){
			itChroms = new LinkedList<>();
			String[] parts = ((String)parameters.getParameterAt(10).getValue()).split(",");
			for(int i=0;i<parts.length;i++){
				itChroms.add(parts[i]);
			}
		}
		
		double perc = (double) parameters.getParameterAt(11).getValue();
		
		
		
		FeatureReader reader = new FeatureReader(numBins, labelsFile, dnaseFile, motifFiles);
		IterativeTraining training = new IterativeTraining(reader, threads, FeatureReader.getSizes(faiFile, bin));
		
		GenDisMixClassifier[] cls = training.iterativeTraining(iterations, trainChroms, itChroms, perc, numBefore, numAfter);
		
		ClassifiersWithInfo info = new ClassifiersWithInfo(cls, numBins, bin, numBefore, numAfter, motifFiles.length);
		
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, info, "classifiers");
		
		File f = File.createTempFile("cls", ".xml");
		f.deleteOnExit();
		
		FileManager.writeFile(f, xml);
		
		TextResult tr = new TextResult("Classifiers", "The trained classifiers", new FileParameter.FileRepresentation(f.getAbsolutePath()), "xml", getToolName(), null, true);
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );
	}

	@Override
	public String getToolName() {
		return "Iterative Training";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "itrain";
	}

	@Override
	public String getDescription() {
		return "performs iterative training on input labels and feature files";
	}

	@Override
	public String getHelpText() {
		return "**Iterative Training** performs an iterative training with the specified number of iterations to obtain a series of classifiers that may be used for predictions in the same cell type or "
				+ "in other cell types based on a corresponding set of feature files. The tool requires as input labels for the training chromosomes, a chromatin accessibility feature file and a set of motif "
				+ "feature files. From the labels, an initial set of training regions is extracted containing all positive examples labeled as *S* (summit) and a sub-sample of negative examples of regions labeled "
				+ "as *U* (unbound). During the iterations, the initial negative examples are complemented with additional negatives obtaining large binding probabilities, i.e., putative false positive predictions. "
				+ "As these additional negative examples are derived from predictions of the current set of classifiers, the number of bins used for aggregation needs to be specified and should be identical to those "
				+ "used for predictions later. Training chromosomes and chromosomes used for predictions in the iterative training may be specified, as well as the percentile of the scores of positive (i.e., summit or bound regions) "
				+ "that should be used to identify putative false positives. The specified bin width must be identical to the bin width specified when computing the corresponding feature files. Feature vectors for training "
				+ "regions may span several adjacent bins as specified by the bin width parameter. Output is an XML file containing the set of trained classifiers. "
				+ "This output file together with a protocol of the tool run is saved to the specified output directory.";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}

}
