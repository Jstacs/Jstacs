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
import de.jstacs.parameters.EnumParameter;
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
import projects.encodedream.UnsupervisedTraining;
import projects.encodedream.UnsupervisedTraining.Init;
import projects.encodedream.UnsupervisedTraining.Select;

public class UnsupervisedTrainingTool implements JstacsTool {

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
		
		
		pars.add(new FileParameter("FAI of genome", "FastA index file of the genome", "fai", true));
		
		try {
			pars.add(new SimpleParameter(DataType.INT, "Bin width", "The width of the genomic bins", true, new NumberValidator<Integer>(1, 1000), 50));
			pars.add(new SimpleParameter(DataType.INT, "Number of bins", "The number of adjacent bins", true, new NumberValidator<Integer>(1, 20), 5));
			pars.add(new SimpleParameter(DataType.INT, "Iterations", "The number of iterations of the interative training", true, new NumberValidator<Integer>(1, 20), 5));
			
			pars.add(new SimpleParameter(DataType.STRING, "Training chromosomes", "Training chromosomes, separated by commas", false));
			
			pars.add(new SimpleParameter(DataType.DOUBLE, "Percentage", "Percentage of positive training examples", true, new NumberValidator<Double>(0.0, 1.0), 0.01));
			pars.add(new SimpleParameter(DataType.DOUBLE, "Factor", "Weight on previous values when computing weights", true, new NumberValidator<Double>(0.0, 100.0), 1.0));
			pars.add(new EnumParameter(UnsupervisedTraining.Init.class, "Initialization of weights", true));
			pars.add(new EnumParameter(UnsupervisedTraining.Select.class, "Selection of next training data set", true));
		} catch (ParameterException e) {
			e.printStackTrace();
		}
		
		
		return new ToolParameterSet(getShortName(),pars.toArray(new Parameter[0]));
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		int numBins = (int) parameters.getParameterAt(4).getValue();
		int bin = (int) parameters.getParameterAt(3).getValue();
		String dnaseFile = (String) parameters.getParameterAt(0).getValue();
		
		int num = ((ExpandableParameterSet) parameters.getParameterAt(1).getValue() ).getNumberOfParameters();
		String[] motifFiles = new String[num];
		for(int i=0;i<num;i++){
			motifFiles[i] = (String) ((ParameterSet)((ExpandableParameterSet) parameters.getParameterAt(1).getValue() ).getParameterAt(i).getValue()).getParameterAt(0).getValue();
		}
		
		String faiFile = (String) parameters.getParameterAt(2).getValue();
		int iterations = (int) parameters.getParameterAt(5).getValue();
		LinkedList<String> trainChroms = null;
		if(((String)parameters.getParameterAt(6).getValue()) != null && ((String)parameters.getParameterAt(6).getValue()).length()>0){
			trainChroms = new LinkedList<>();
			String[] parts = ((String)parameters.getParameterAt(6).getValue()).split(",");
			for(int i=0;i<parts.length;i++){
				trainChroms.add(parts[i]);
			}
		}
		
		double frac = (double) parameters.getParameterAt(7).getValue();
		double cons = (double) parameters.getParameterAt(8).getValue();
		Init init = (Init) parameters.getParameterAt(9).getValue();
		Select select = (Select) parameters.getParameterAt(10).getValue();
		
		
		FeatureReader reader = new FeatureReader(numBins, null, dnaseFile, motifFiles);
		
		
		UnsupervisedTraining training = new UnsupervisedTraining(reader, threads, FeatureReader.getSizes(faiFile, bin), init, select);
		
		
		GenDisMixClassifier[] cls = training.iterativeTraining(iterations, trainChroms, frac, cons);
		
		ClassifiersWithInfo info = new ClassifiersWithInfo(cls, numBins, bin, 1, numBins-1, motifFiles.length);
		
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
		return "Unsupervised Training";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "unsuper";
	}

	@Override
	public String getDescription() {
		return "performs unsupervised training on feature files";
	}

	@Override
	public String getHelpText() {
		return "";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}

}
