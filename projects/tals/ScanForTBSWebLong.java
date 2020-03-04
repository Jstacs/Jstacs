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

import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.util.Arrays;
import java.util.Comparator;

import projects.tals.TBSScanner.ResultList;
import projects.tals.TBSScannerLong.TBSScannerLongParameterSet;
import de.jstacs.DataType;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotationParser;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.FileManager;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.DataSetResult;
import de.jstacs.results.ImageResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.DifferentiableStatisticalModelWrapperTrainSM;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor.Protocol;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.PFMComparator;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SeqLogoPlotter;

/**
 * Main class of TALgetterLongWeb.jar
 * 
 * @author Jan Grau
 *
 */
public class ScanForTBSWebLong {

	
	/**
	 * Main method for creating Galaxy tool xml and for running the Galaxy tool.
	 * 
	 * @param args for creation of Galaxy tool xml: --create <path_to.xml>, otherwise parameters provided by Galaxy
	 * @throws Exception if something went wrong
	 */
	public static void main( String[] args ) throws Exception {
				
		
		
		TBSScannerLongParameterSet params = new TBSScannerLongParameterSet();
		SelectionParameter mod2 = new SelectionParameter( DataType.INT, new String[]{"TALgetter", "TALgetter13"}, new Integer[]{1,2}, "Model type", 
				"TALgetter is the default model that uses individual binding specificities for each RVD. " +
				"TALgetter13 uses binding specificities that only depend on amino acid 13, i.e., the second amino acid of the repat." +
				"While TALgetter is recommended in most cases, the use of TALgetter13 may be beneficial if you search for target sites of TAL effector with many rare RVDs, for instance YG, HH, or S*.", true );
		mod2.setDefault( "TALgetter" );
		
		SelectionParameter train = new SelectionParameter( DataType.PARAMETERSET, new String[]{"Use standard model", "Use previously trained model", "Train model on training data"}, new ParameterSet[]{
					                                                                                                                                         new SimpleParameterSet( mod2 ),
					                                                                                                                                         new SimpleParameterSet( new FileParameter( "Model", "Choose a TALgetter model from your history", "xml", true ) ),
					                                                                                                                                         new SimpleParameterSet(mod2, new FileParameter( "Training data", "The training data, annotated FastA format. The required format is described in the help section, where we also provide the data set used to train the default models.", "fasta", true ) )
					}, "Model training", "You can either use the standard TALgetter model, re-use a TALgetter model that has already been trained on given pairs of TAL effectors and target sites, or you can provide your own training data. ", true );
		
		
		params.addParameter( params.getNumberOfParameters(), train );
		
		boolean[] line = new boolean[]{true,true,true,false,false,true,true};
		
		SimpleParameterSet ps = new SimpleParameterSet( params.getAllParameters() );
		
		GalaxyAdaptor ga = new GalaxyAdaptor( ps,null, line,"TALgetterLong", "TALgetterLong (TAL effector target site finder) is a variant of TALgetter that is specifically designed for large input data sets.", "1.0", "java -Xms256M -Xmx2G -jar "+System.getProperty( "user.dir" )+System.getProperty( "file.separator" )+"TALgetterWebLong.jar", "jobname", true );
		ga.setHelp( FileManager.readInputStream( ScanForTBSWebLong.class.getClassLoader().getResourceAsStream( "projects/tals/helplong.txt" ) ).toString() );
		
		if(!ga.parse( args, false )){
			System.exit( 1 );
		}	
		
		AbstractDifferentiableStatisticalModel model = null;
		
		if(((SelectionParameter)params.getParameterForName( "Model training" )).getSelected() != 1){
			mod2 = (SelectionParameter)((SimpleParameterSet)params.getParameterForName( "Model training" ).getValue()).getParameterAt( 0 );

			if(mod2.getValue().equals( 1 )){
				model = (AbstractDifferentiableStatisticalModel)XMLParser.extractObjectForTags( FileManager.readInputStream( ScanForTBSWebLong.class.getClassLoader().getResourceAsStream( "projects/tals/talfinder_obg2_hyp_bg.xml" ) ), "model" );
			}else{
				model = (AbstractDifferentiableStatisticalModel)XMLParser.extractObjectForTags( FileManager.readInputStream( ScanForTBSWebLong.class.getClassLoader().getResourceAsStream( "projects/tals/talfinder_obg2_hyp_bg_map.xml" ) ), "model" );
			}
		}
		
		Protocol prot = ga.getProtocol( false );
		ByteArrayOutputStream baos = prot.getOutputStream();

		String[] alph={"NI","NG","NN","NS","N*","ND","NK","NC","NV","NA","NH","HD","HG","HA","H*","HH","HI","HN","S*","SN","SS","IG","YG","NP","NT","IS"};
		AlphabetContainer alphabetsRVD= new AlphabetContainer( new DiscreteAlphabet(true,alph) );
		
		if(((SelectionParameter)params.getParameterForName( "Model training" )).getSelected() == 1){
			FileParameter fp = (FileParameter)((SimpleParameterSet)params.getParameterForName( "Model training" ).getValue()).getParameterAt( 0 );
			
			String fn = (String)fp.getValue();
			StringBuffer sb = FileManager.readFile( fn );
			model = new TALgetterDiffSM( sb );
		}else if(((SelectionParameter)params.getParameterForName( "Model training" )).getSelected() == 2){
			
			
			String filename = (String) ((SimpleParameterSet)params.getParameterForName( "Model training" ).getValue()).getParameterAt( 1 ).getValue();
			
			DataSet trainDs = new DNADataSet( filename, '>', new ReferenceSequenceAnnotationParser( "seq", alphabetsRVD, ":", ";", "-" ) );
			
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
			
			model =  (AbstractDifferentiableStatisticalModel)trainer.getFunction();
			
			
			String[] alph2 = alph.clone();
			Arrays.sort( alph2, new Comparator<String>() {@Override
			public int compare( String o1, String o2 ) {
				int c = Character.valueOf( o1.charAt( 1 ) ).compareTo( Character.valueOf( o2.charAt( 1 ) ) );
				if(c == 0){
					return Character.valueOf( o1.charAt( 0 ) ).compareTo( Character.valueOf( o2.charAt( 0 ) ) );
				}else{
					if(o1.charAt( 1 ) == '*'){
						return 1;
					}else if(o2.charAt( 1 ) == '*'){
						return -1;
					}else{
						return c;
					}
				}
			}} );
			StringBuffer modSeq = new StringBuffer();
			for(int i=0;i<alph2.length;i++){
				modSeq.append( alph2[i] );
				if(i < alph2.length-1){
					modSeq.append( "-" );
				}
			}
			
			Pair<double[][],double[]> specNImpAll = ((TALgetterDiffSM)model).getSpecificitiesAndImportances( Sequence.create( alphabetsRVD, modSeq.toString(), "-" ) );
			
			int height = SeqLogoPlotter.getHeight( 750, specNImpAll.getFirstElement() );
			
			BufferedImage img = SeqLogoPlotter.plotTALgetterLogoToBufferedImage( height, specNImpAll.getFirstElement(), specNImpAll.getSecondElement(), ("0-"+modSeq.toString()).split( "-" ) );
			
			ga.addResult( new ImageResult( "Model logo", "Logo plot representing specificities and importances learned from the training data", img ), false, true );

			
			prot.appendHeading( "Training TALgetter model" );
			prot.append( "Trained on "+weights.length+" pairs of TAL effector and target site.<br /><br />" );
			prot.append( "Model parameters: <br />" );
			prot.append( model.toString().replaceAll( "\\n", "<br />" ) );
			prot.append( "<br /><br />" );
			ga.addResult( new StorableResult( "TALgetter model", "TALgetter model", model ), true, false );
		}
		
		
		((TALgetterDiffSM)model).fix();//TODO
		
		
		
		BufferedReader dsRead = params.getInputSequences();
		
		prot.appendHeading( "Search" );
		prot.append( "Searching for targets of the TAL effector with RVD sequence<br />" );
		prot.append( "    "+params.getTALSequence()+".<br />" +
				"using "+(model instanceof TALgetter13DiffSM ? "TALgetter13" : "TALgetter")+".<br /><br />" );
		prot.append( "Reporting best "+params.getN()+" target sites.<br /><br />" );
		
		
		ResultList[] rls = TBSScannerLong.scan( (TALgetterDiffSM)model, params, dsRead );
		
		ga.addResult( new ListResult( "Description of output columns", "The output can also be downloaded as a tab-separated file.", null, 
				new ResultSet[]{
				                new ResultSet( new Result[]{new CategoricalResult("Column","","ID"),
				                                            new CategoricalResult( "Description", "", "The ID of the input sequence as given in FastA header" )} ),
	                            new ResultSet( new Result[]{new CategoricalResult("Column","","Position"),
	                                                        new CategoricalResult( "Description", "", "Position in the given sequence" )} ),
				                new ResultSet( new Result[]{new CategoricalResult("Column","","Distance to end"),
				                                            new CategoricalResult( "Description", "", "The distance to the right end of the sequence" )} ),
				                new ResultSet( new Result[]{new CategoricalResult("Column","","Sequence"),
				            				                new CategoricalResult( "Description", "", "The sequence of the target site site" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Matches"),
				            				                new CategoricalResult( "Description", "", "Categories of matches per position: M/m first position match/mismatch; | match; : weak match; x mismatch" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Score"),
				            				                new CategoricalResult( "Description", "", "Score returned by the model" )} )
				            				                                            
				            				                                            
				            				                                            
				            				                                            
				}
		), false, true );
		

		DataSet bs = rls[0].getBindingSites();
		
		
		
		ResultSet[] res = rls[0].toArray( );
		
		if(res.length > 0){
			ga.addResult( new ListResult( "Predictions - tabular", "", null, res ), true, false );
		}
		ga.addResult( new ListResult( "Predictions", "The predictions are also available as a tab-seperated file and as a FastA-file containing all predicted target sites from your history.", null, rls[1].toArray( ) ), false, true );
		
		
		if(res.length > 0){
			ga.addResult( new DataSetResult( "Binding sites - FastA", "The top "+res.length+" binding sites" , bs ), true, false );
		}
				
		
		Pair<double[][],double[]> specNImp = ((TALgetterDiffSM)model).getSpecificitiesAndImportances( Sequence.create( alphabetsRVD, params.getTALSequence(), "-" ) );
		
		int height = SeqLogoPlotter.getHeight( 750, specNImp.getFirstElement() );
		
		BufferedImage img = SeqLogoPlotter.plotTALgetterLogoToBufferedImage( height, specNImp.getFirstElement(), specNImp.getSecondElement(), ("0-"+params.getTALSequence()).split( "-" ) );
		
		ga.addResult( new ImageResult( "Theoretical target site logo", "Logo plot representing specificities and importances for the input TAL effector according to model parameters", img ), false, true );
		
		if(res.length > 0){
			double[][] pfm = PFMComparator.getPFM( bs );
			for(int i=0;i<pfm.length;i++){
				Normalisation.sumNormalisation( pfm[i] );
			}

			BufferedImage bsLogo = SeqLogoPlotter.plotLogoToBufferedImage( height, pfm );

			ga.addResult( new ImageResult( "Predicted target site logo", "Sequence logo of the predicted target sites, depends on the threshold on the p-values and the maximum number of reported target sites", bsLogo ), false, true );
		}
		prot.appendHeading( "Finished...<br /><br />" );
		
		prot.appendHeading( "Result" );
		prot.append( "Found "+rls[0].getNumberOfResults()+" target sites with scores between "+rls[0].getBestScore()+" and "+rls[0].getWorstScore()+"." );
		
		ga.writeOutput();
	}
	
	
}
