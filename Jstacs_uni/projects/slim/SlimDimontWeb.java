package projects.slim;

import java.io.ByteArrayOutputStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Locale;

import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.DataSet;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor.Protocol;
import de.jstacs.utils.SafeOutputStream;
import projects.dimont.ThresholdedStrandChIPper;


public class SlimDimontWeb extends SlimDimont {

	/**
	 * @param args
	 */
	public static void main( String[] args ) throws Exception {
		
		
		SlimDimontWebParameterSet params = new SlimDimontWebParameterSet();
		
		boolean[] lines = new boolean[params.getNumberOfParameters()];
		lines[0] = lines[5] = lines[6] = lines[9] = true;
		
		
		GalaxyAdaptor ga = new GalaxyAdaptor( params,lines,"SlimDimont", "- a universal tool for detecting dependencies from ChIP and PBM data.", "0.1", "java -Xms256M -Xmx2G -jar "+System.getProperty( "user.dir" )+System.getProperty( "file.separator" )+"SlimDimontWeb.jar", "jobname" );
		ga.setHelp( FileManager.readInputStream( SlimDimontWeb.class.getClassLoader().getResourceAsStream( "projects/slim/help.txt" ) ).toString() );
		
		if(!ga.parse( args, true )){
			System.exit( 1 );
		}	
		
		int threads = ga.getThreads();
		
		DataSet data = params.getInputSequences();
		int motifLength = (Integer)params.getParameterForName( "Motif width" ).getValue();
		int restarts = (Integer)params.getParameterForName( "Starts" ).getValue();
		int fgOrder = 0;
		
		int selected = ((SelectionParameter)params.getParameterForName( "Model type" )).getSelected();
		if(selected == 0){
			fgOrder = (Integer)((ParameterSet)((SelectionParameter)params.getParameterForName( "Model type" )).getValue()).getParameterAt( 0 ).getValue();
		}else if(selected == 1){
			fgOrder = -motifLength;
		}else{
			fgOrder = -(Integer)((ParameterSet)((SelectionParameter)params.getParameterForName( "Model type" )).getValue()).getParameterAt( 0 ).getValue();
		}
		
		int bgOrder = (Integer)params.getParameterForName( "Markov order of background model" ).getValue();
		String position = (String)params.getParameterForName( "Position tag" ).getValue();
		String value = (String)params.getParameterForName( "Value tag" ).getValue();
		double sd = (Double)params.getParameterForName( "Standard deviation" ).getValue();
		String weightingFactor = (String)params.getParameterForName( "Weighting factor" ).getValue();
		double ess = (Double)params.getParameterForName( "Equivalent sample size" ).getValue();
		boolean delete = (Boolean)params.getParameterForName( "Delete BSs from profile" ).getValue();
		boolean modify = (Boolean)params.getParameterForName( "Adjust for shifts" ).getValue();
		
		Protocol prot = ga.getProtocol( false );
		ByteArrayOutputStream baos = prot.getOutputStream();
		
		
		Result[][] res = run( data,motifLength,restarts,fgOrder,bgOrder,position,value,weightingFactor,ess,delete,threads, SafeOutputStream.getSafeOutputStream( baos ),sd,modify);

		
		
		
		NumberFormat nf = DecimalFormat.getInstance( Locale.ENGLISH );
		nf.setMaximumFractionDigits( 3 );
		nf.setMinimumFractionDigits( 3 );
		
		for(int i=0;i<res.length;i++){
			prot.append( "\n++++++++++++++++++++++++++++++++++++++++++++++\n\n" );
			prot.appendHeading( "Motif model "+(i+1) );
			prot.append( ((ThresholdedStrandChIPper)((GenDisMixClassifier)((StorableResult)res[i][0]).getResultInstance()).getDifferentiableSequenceScore( 0 ) ).toHtml( nf ) );
			ga.addResult( res[i][0], true, false );
			if(res[i].length > 2){
				ga.addResult( res[i][2], false, true);
				ga.addResult( res[i][3], false, true);
				ga.addResult( res[i][4], false, true);
			}
		}
		
		for(int i=0;i<res.length;i++){
			ga.addResult( res[i][1], true, false );
		}
		
		
		ga.addResult( new ListResult( "Description of columns of binding site output (see history)", "You can download the predictions for the motifs discovered as tab-separated file from the history.", null, 
				new ResultSet[]{
				                new ResultSet( new Result[]{new CategoricalResult("Column","","Sequence index"),
				                                            new CategoricalResult( "Description", "", "The index of the sequence" )} ),
	                            new ResultSet( new Result[]{new CategoricalResult("Column","","Position"),
	                                                        new CategoricalResult( "Description", "", "The start position of predicted binding site (BS) within the sequence" )} ),
				                new ResultSet( new Result[]{new CategoricalResult("Column","","Strand"),
				                                            new CategoricalResult( "Description", "", "The strand of the predicted BS" )} ),
				                new ResultSet( new Result[]{new CategoricalResult("Column","","p-value"),
				            				                new CategoricalResult( "Description", "", "The p-value of the predicted BS" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","-log10(p-value)"),
				            	                            new CategoricalResult( "Description", "", "The negative logarithm of the p-value of the predicted BS" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Score"),
				            				            	new CategoricalResult( "Description", "", "The model score of the predicted BS" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Binding site"),
				            				                new CategoricalResult( "Description", "", "The binding site as in the sequence" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Adjusted binding site"),
				            				                new CategoricalResult( "Description", "", "The binding site in predicted orientation" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Signal"),
				            				            	new CategoricalResult( "Description", "", "The signal of the sequence annotation" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Sequence annotation"),
				            				                new CategoricalResult( "Description", "", "The annotation of the original sequence" )} )
			}
		), false, true );
		
		ga.writeOutput();
		
		
	}

}
