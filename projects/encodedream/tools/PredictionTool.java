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
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;

import de.jstacs.DataType;
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
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import projects.encodedream.ClassifiersWithInfo;
import projects.encodedream.FeatureReader;
import projects.encodedream.Predictor;

public class PredictionTool implements JstacsTool {

	@Override
	public ToolParameterSet getToolParameters() {
		
		LinkedList<Parameter> pars = new LinkedList<>();
		
		pars.add(new FileParameter("Classifiers", "The classifiers trained by iterative training", "xml", true));
		
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
			pars.add(new SimpleParameter(DataType.STRING, "Prediction chromosomes", "Prediction chromosomes, separated by commas", false));
			
			pars.add(new SimpleParameter(DataType.INT,"Aggregation: bins before","Number of bins before the current one considered for aggregation.",false));
			pars.add(new SimpleParameter(DataType.INT,"Aggregation: bins after","Number of bins after the current one considered for aggregation.",false));
			
			pars.add(new SimpleParameter(DataType.INT,"Number of classifiers","Use only the first k (last k for negative values) classifiers for predictions.",false));
			
		} catch (ParameterException e) {
			e.printStackTrace();
		}
		
		return new ToolParameterSet(getShortName(),pars.toArray(new Parameter[0]));
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		ClassifiersWithInfo cls = (ClassifiersWithInfo) XMLParser.extractObjectForTags(FileManager.readFile((String)parameters.getParameterAt(0).getValue()), "classifiers");
		
		String dnaseFile = (String) parameters.getParameterAt(1).getValue();
		
		int num = ((ExpandableParameterSet) parameters.getParameterAt(2).getValue() ).getNumberOfParameters();
		String[] motifFiles = new String[num];
		for(int i=0;i<num;i++){
			motifFiles[i] = (String) ((ParameterSet)((ExpandableParameterSet) parameters.getParameterAt(2).getValue() ).getParameterAt(i).getValue()).getParameterAt(0).getValue();
		}
		
		if(motifFiles.length != cls.getNumberOfMotifs()){
			throw new Exception("Not the same number of motifs");
		}
		
		String faiFile = (String) parameters.getParameterAt(3).getValue();
		
		HashMap<String,Integer> sizes = FeatureReader.getSizes(faiFile, cls.getBinWidth());
		
		LinkedList<String> predChroms = new LinkedList<String>(sizes.keySet());
		Collections.sort(predChroms);
		if(((String)parameters.getParameterAt(4).getValue()) != null && ((String)parameters.getParameterAt(4).getValue()).length()>0){
			predChroms = new LinkedList<>();
			String[] parts = ((String)parameters.getParameterAt(4).getValue()).split(",");
			for(int i=0;i<parts.length;i++){
				predChroms.add(parts[i]);
			}
		}
		
		Integer binsBefore = (Integer) parameters.getParameterAt(5).getValue();
		if(binsBefore == null){
			binsBefore = cls.getBinsBefore();
		}
		Integer binsAfter = (Integer) parameters.getParameterAt(6).getValue();
		if(binsAfter == null){
			binsAfter = cls.getBinsAfter();
		}
		Integer numClass = (Integer) parameters.getParameterAt(7).getValue();
		if(numClass != null){
			cls.limitClassifiers(numClass);
		}
		System.out.println(binsBefore+" "+binsAfter);
		
		
		FeatureReader reader = new FeatureReader(cls.getNumBins(), null, dnaseFile, motifFiles);
		Predictor pred = new Predictor(cls.getClassifiers(),reader,binsBefore,binsAfter);
		
		File f = pred.predict(sizes,predChroms);
		
		TextResult tr = new TextResult("Predictions", "Predictions of binding probabilities in tabular format", new FileParameter.FileRepresentation(f.getAbsolutePath()), "tsv.gz", getToolName(), null, true);
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );
		
	}

	@Override
	public String getToolName() {
		return "Prediction";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "predict";
	}

	@Override
	public String getDescription() {
		return "predicts binding probabilities of genomic regions";
	}

	@Override
	public String getHelpText() {
		return "**Prediction** predicts binding probabilities of genomic regions as specified during training of the set of classifiers in iterative training. As input, Prediction requires a set of trained classifiers "
				+ "in XML format, the same (type of) feature files as used in training (motif files must be specified in the same order!). In addition, the chromosomes for which predictions are made may be specified, "
				+ "and the number of bins used for aggregation may be specified to deviate from those used during training. If these bin numbers are not specified, those from the training run are used. Finally, it is possible "
				+ "to restrict the number of classifiers considered to the first n ones. Output is provided as a gzipped file *Predictions.tsv.gz* with columns chromosome, start position, binding probability. "
				+ "This output file together with a protocol of the tool run is saved to the specified output directory.";
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
