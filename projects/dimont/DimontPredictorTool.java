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

package projects.dimont;

import java.io.IOException;
import java.io.StringReader;
import java.util.Date;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.SparseSequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.SplitSequenceAnnotationParser;
import de.jstacs.io.FileManager;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.motifDiscovery.SignificantMotifOccurrencesFinder;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.ImageResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.JstacsTool.ResultEntry;
import de.jstacs.utils.PFMComparator;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SeqLogoPlotter;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.SeqLogoPlotter.SeqLogoPlotGenerator;
import projects.motifComp.FindPWMsAndClusters;

public class DimontPredictorTool implements JstacsTool {

	public DimontPredictorTool() {
	}

	@Override
	public ParameterSet getToolParameters() {
		
		LinkedList<Parameter> parameters = new LinkedList<Parameter>();
		
		parameters.add(new FileParameter("Input file", "The file name of the file containing the input sequences in annotated FastA format (see readme)", "fasta,fa,fas", true));
		
		parameters.add(new FileParameter("Dimont classifier", "The classifier from the Dimont output for one motif", "xml", true));
		
		try{
		parameters.add( new SimpleParameter( DataType.STRING, "Value tag", "The tag for the value information in the FastA-annotation of the input file", true, "signal" ) );
		
		parameters.add( new SimpleParameter( DataType.STRING, "Weighting factor", "The value for weighting the data; either a value between 0 and 1, or a description relative to the standard deviation (e.g. +4sd)", true, "" + 0.2 ) );
		
		parameters.add( new SimpleParameter( DataType.DOUBLE, "p-value", "The maximum p-value allowed for predicted binding sites", true, new NumberValidator<Double>( 0.0, 1.0 ), 1E-3 ) );
		
		}catch(Exception e){
			throw new RuntimeException();
		}
		
		return new SimpleParameterSet(parameters.toArray(new Parameter[0]));
		
	}

	@Override
	public ToolResult run(ParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		
		
		DataSet data = SparseSequence.getDataSet( DNAAlphabetContainer.SINGLETON, new SparseStringExtractor( 
				new StringReader(
						((FileParameter)parameters.getParameterAt(0)).getFileContents().getContent())
				, '>', "", new SplitSequenceAnnotationParser(":",";") ) );
		
		
		GenDisMixClassifier cl = new GenDisMixClassifier( new StringBuffer(((FileParameter)parameters.getParameterAt(1)).getFileContents().getContent()) );
		
		ThresholdedStrandChIPper model = (ThresholdedStrandChIPper)cl.getDifferentiableSequenceScore( 0 );
		
		String value = parameters.getParameterAt(2).getValue().toString();
		
		String weightingFactor = parameters.getParameterAt(3).getValue().toString();
		
		double pval = (Double) parameters.getParameterAt(4).getValue();
		
		double[][] weights = new double[2][data.getNumberOfElements()];
		
		double[] raw = weights[0].clone();
		
		//read annotation
		for( int j = 0; j < weights[0].length; j++ ) {
			Sequence seq = data.getElementAt(j);
			SequenceAnnotation[] seqAn = seq.getAnnotation();
			for( int i = 0; i < seqAn.length; i++ ) {
				if( seqAn[i].getType().equals(value) ) {
					raw[j] = Double.parseDouble( seqAn[i].getIdentifier() );
				}
			}
		}
		
		//create weights
		double wf;
		if( weightingFactor.endsWith("sd") ) {
			double h = Double.parseDouble( weightingFactor.substring(0,weightingFactor.length()-2) );
			double meanRaw = ToolBox.sum(raw) / raw.length;
			double sdRaw = 0;
			for( int i = 0; i < raw.length; i++ ) {
				sdRaw += (raw[i]-meanRaw) * (raw[i]-meanRaw);
			}
			sdRaw = Math.sqrt( sdRaw/raw.length );
			h = meanRaw + h*sdRaw;
			double anz = 0;
			for( int i = 0; i < raw.length; i++ ) {
				if( raw[i] >= h ) {
					anz++;
				}
			}
			anz=Math.max(50,anz);
			wf = anz/raw.length;
		} else {
			wf = Double.parseDouble( weightingFactor );
		}

		weights[0] = Interpolation.getWeight( data, raw, wf, Interpolation.RANK_LOG );
		weights[1] = Interpolation.getBgWeight( weights[0] );
		
	
		
		SignificantMotifOccurrencesFinder smof = new SignificantMotifOccurrencesFinder( model, data, weights[1], pval );
		
		Pair<double[][][],int[][]> pair = smof.getPWMAndPositions( 0, data, weights[0], 0, 0 );
		
		LinkedList<Result> result = new LinkedList<Result>();
		
		
		
		result.add(DimontTool.getListResult(data, weights[0],pair, model.getMotifLength( 0 ), 0 ));
		
		double[][] pwm = pair.getFirstElement()[0];
		
		if(!Double.isNaN( pwm[0][0] )){
			
			try{
				//int height = SeqLogoPlotter.getHeight( 750, pwm );
				
				result.add(new PlotGeneratorResult( "Sequence logo", "Sequence logo of the motif ", 
						new SeqLogoPlotGenerator(pwm, 1000), false));
				
				result.add(new PlotGeneratorResult( "Sequence logo (rc)", "Sequence logo of the reverse complement of the motif", 
						new SeqLogoPlotGenerator(PFMComparator.getReverseComplement( DNAAlphabet.SINGLETON, pwm ), 1000), false));
				
									}catch(Exception e){
				
			}catch(InternalError err){
				
			}
			
			
			
		}
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(result), parameters, getToolName(), new Date(System.currentTimeMillis()));
		
		
	}

	@Override
	public String getToolName() {
		return "Dimont Predictor";
	}

	@Override
	public String getToolVersion() {
		return "1.2";
	}

	@Override
	public String getShortName() {
		return "predict";
	}

	@Override
	public String getDescription() {
		return "for predicting binding sites using a Dimont model";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( FindPWMsAndClusters.class.getClassLoader().getResourceAsStream( "projects/dimont/helpPredictor.txt" ) ).toString();
		} catch ( IOException e ) {
			e.printStackTrace();
			return "";
		}
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[]{
				new ResultEntry(ListResult.class, null, "Predictions for motif 1")
		};
	}

}
