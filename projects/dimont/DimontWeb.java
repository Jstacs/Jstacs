package projects.dimont;

import java.io.ByteArrayOutputStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Locale;

import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.DataSet;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.Parameter;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;
import de.jstacs.utils.SafeOutputStream;
import de.jstacs.utils.galaxy.GalaxyAdaptor;
import de.jstacs.utils.galaxy.GalaxyAdaptor.Protocol;


public class DimontWeb {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main( String[] args ) throws Exception {
		
		DimontWebParameterSet params = new DimontWebParameterSet();
		
		boolean[] lines = new boolean[params.getNumberOfParameters()];
		lines[0] = lines[5] = lines[6] = lines[9] = true;
		
		
		GalaxyAdaptor ga = new GalaxyAdaptor( params,lines,"Dimont", "- a universal tool for de-novo motif discovery (beta).", "0.1", "java -Xms256M -Xmx2G -jar "+System.getProperty( "user.dir" )+System.getProperty( "file.separator" )+"DimontWeb.jar", "jobname" );
		ga.setHelp( FileManager.readInputStream( DimontWeb.class.getClassLoader().getResourceAsStream( "projects/dimont/help.txt" ) ).toString() );
		
		if(!ga.parse( args )){
			System.exit( 1 );
		}	
		
		DataSet data = params.getInputSequences();
		int motifLength = (Integer)params.getParameterForName( "Initial motif width" ).getValue();
		int restarts = (Integer)params.getParameterForName( "Starts" ).getValue();
		int fgOrder = (Integer)params.getParameterForName( "Markov order of motif model" ).getValue();
		int bgOrder = (Integer)params.getParameterForName( "Markov order of background model" ).getValue();
		String position = (String)params.getParameterForName( "Position tag" ).getValue();
		String value = (String)params.getParameterForName( "Value tag" ).getValue();
		double sd = (Double)params.getParameterForName( "Standard deviation" ).getValue();
		String weightingFactor = (String)params.getParameterForName( "Weighting factor" ).getValue();
		double ess = (Double)params.getParameterForName( "Equivalent sample size" ).getValue();
		boolean delete = (Boolean)params.getParameterForName( "Delete BSs from profile" ).getValue();
		int threads = 1;
		
		Protocol prot = ga.getProtocol( false );
		ByteArrayOutputStream baos = prot.getOutputStream();
		
		Result[][] res = Dimont.run( data,motifLength,restarts,fgOrder,bgOrder,position,value,sd,weightingFactor,ess,delete,threads, SafeOutputStream.getSafeOutputStream( baos ));
		
		NumberFormat nf = DecimalFormat.getInstance( Locale.ENGLISH );
		nf.setMaximumFractionDigits( 3 );
		nf.setMinimumFractionDigits( 3 );
		
		for(int i=0;i<res.length;i++){
			prot.append( "\n++++++++++++++++++++++++++++++++++++++++++++++\n\n" );
			prot.appendHeading( "Motif model "+(i+1) );
			prot.append( ((ThresholdedStrandChIPper)((GenDisMixClassifier)((StorableResult)res[i][0]).getResultInstance()).getDifferentiableSequenceScore( 0 ) ).toHtml( nf ) );
			ga.addResult( res[i][0], true, false );
			if(res[i].length > 3){
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
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Binding site"),
				            				                new CategoricalResult( "Description", "", "The binding site as in the sequence" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Adjusted binding site"),
				            				                new CategoricalResult( "Description", "", "The binding site in predicted orientation" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Sequence annotation"),
				            				                new CategoricalResult( "Description", "", "The annotation of the original sequence" )} )
			}
		), false, true );
		
		ga.writeOutput();
		
	}

}
