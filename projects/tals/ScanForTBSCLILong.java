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
package projects.tals;

import java.io.BufferedReader;
import java.util.Arrays;

import projects.tals.TBSScanner.ResultList;
import projects.tals.TBSScannerLong.TBSScannerLongParameterSet;
import de.jstacs.DataType;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotationParser;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.FileManager;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.ParameterSetTagger;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.ListResult;
import de.jstacs.results.ResultSet;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.DifferentiableStatisticalModelWrapperTrainSM;

/**
 * Main class of TALgetterLong.jar, a variant designed for large input data sets.
 * 
 * @author Jan Grau
 *
 */
public class ScanForTBSCLILong {
	
	
	/**
	 * Program for scanning input sequences for TAL effector target sites on the command line.
	 * 
	 * @param args parameters in key=value format
	 * @throws Exception if something went wrong
	 */
	public static void main( String[] args ) throws Exception {
		
		TBSScannerLongParameterSet params = new TBSScannerLongParameterSet();
		
		SelectionParameter mod = new SelectionParameter( DataType.STRING, new String[]{"TALgetter", "TALgetter13"}, new String[]{"TALgetter", "TALgetter13"}, "Model type", 
				"TALgetter is the default model that uses individual binding specificities for each RVD. " +
				"TALgetter13 uses binding specificities that only depend on amino acid 13, i.e., the second amino acid of the repat." +
				"While TALgetter is recommended in most cases, the use of TALgetter13 may be beneficial if you search for target sites of TAL effector with many rare RVDs, for instance YG, HH, or S*.", true );
		mod.setDefault( "TALgetter" );
		
		SimpleParameterSet tempPars = new SimpleParameterSet(
				new SimpleParameter( DataType.STRING, "Input sequences", "The sequences to scan for TAL effector target sites, FastA", true ),
				params.getParameterForName( "RVD sequence" ),
				mod,
				params.getParameterForName( "Maximum number of target sites" ),
				new SimpleParameter(DataType.STRING, "Training data", "The input data to use for training the model, annotated FastA", false)
				);
		
		ParameterSetTagger pst = new ParameterSetTagger( 
				new String[]{"input","rvd","model","top","train"}, 
				tempPars);
		
		pst.fillParameters( "=", args );
		
		
		System.err.println( "parameters:" );
		System.err.println( pst );
		System.err.println("_________________________________");
		if( !tempPars.hasDefaultOrIsSet() ) {
			System.err.println( "Some of the required parameters are not specified." );
			System.exit( 1 );
		}
		params.setInputPath( (String)tempPars.getParameterAt( 0 ).getValue() );
		
		DifferentiableStatisticalModel model = null;
		
		if(mod.getValue().equals( "TALgetter" )){
			model = (DifferentiableStatisticalModel)XMLParser.extractObjectForTags( FileManager.readInputStream( ScanForTBSCLILong.class.getClassLoader().getResourceAsStream( "projects/tals/talfinder_obg2_hyp_bg.xml" ) ), "model" );		
		}else{
			model = (DifferentiableStatisticalModel)XMLParser.extractObjectForTags( FileManager.readInputStream( ScanForTBSCLILong.class.getClassLoader().getResourceAsStream( "projects/tals/talfinder_obg2_hyp_bg_map.xml" ) ), "model" );
		}
		
		if(tempPars.getParameterForName( "Training data" ).getValue() != null){
			String[] alph={"NI","NG","NN","NS","N*","ND","NK","NC","NV","NA","NH","HD","HG","HA","H*","HH","HI","HN","S*","SN","SS","IG","YG","NP","NT","IS"};
			AlphabetContainer alphabetsRVD= new AlphabetContainer( new DiscreteAlphabet(true,alph) );
			
			DNADataSet trainDs = new DNADataSet( (String)tempPars.getParameterForName( "Training data" ).getValue(), '>', new ReferenceSequenceAnnotationParser( "seq", alphabetsRVD, ":", ";", "-" ) );
			double[] weights = new double[trainDs.getNumberOfElements()];
			Arrays.fill( weights, 1.0 );
			for(int i=0;i<weights.length;i++) {
				Sequence seq = trainDs.getElementAt( i );
				SequenceAnnotation ann = seq.getSequenceAnnotationByType( "weight", 0 );
				if(ann != null){
					weights[i] = Double.parseDouble( ann.getIdentifier() );
				}
			}
			
			DifferentiableStatisticalModelWrapperTrainSM trainer =  new DifferentiableStatisticalModelWrapperTrainSM( model, 1, Optimizer.QUASI_NEWTON_BFGS, new SmallDifferenceOfFunctionEvaluationsCondition( 1E-12 ), 1E-12, 1E-4 );
			trainer.setOutputStream( null );
			trainer.train(trainDs,weights);
			
			model =  trainer.getFunction();
			System.err.println("Estimated parameters: ");
			System.err.println(model);
			System.err.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
		}
		
		
		((TALgetterDiffSM)model).fix();
		
		BufferedReader ds = params.getInputSequences();
		
		System.err.println( "Searching for targets of TAL effector sequence" );
		System.err.println( "\t"+params.getTALSequence()+".\n" );
		System.err.println( "Reporting best "+params.getN()+" target sites.\n\n" );
		
		
		ResultList[] rls = TBSScannerLong.scan( (TALgetterDiffSM)model, params, ds );
		
		ResultSet[] res = rls[0].toArray( );
		
		System.out.println( new ListResult( "Predictions", "", null, res ) );
		
		System.err.println( "Finished..." );
		
	}
	
	
}
