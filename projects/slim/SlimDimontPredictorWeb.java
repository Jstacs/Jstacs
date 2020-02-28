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

package projects.slim;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Locale;

import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.DataSet;
import de.jstacs.io.FileManager;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor.Protocol;
import projects.dimont.ThresholdedStrandChIPper;


public class SlimDimontPredictorWeb {

	/**
	 * @param args
	 */
	public static void main( String[] args ) throws Exception {
		
		
		SlimDimontPredictorWebParameterSet params = new SlimDimontPredictorWebParameterSet();
		
		boolean[] lines = new boolean[params.getNumberOfParameters()];
		lines[0] = lines[1] = lines[4] = true;
		
		
		GalaxyAdaptor ga = new GalaxyAdaptor( params,null,lines,"SlimDimontPredictor", "for predicting binding sites using a Dimont model", "0.1", "java -Xms256M -Xmx2G -jar "+System.getProperty( "user.dir" )+System.getProperty( "file.separator" )+"SlimDimontPredictorWeb.jar", "jobname", true );
		ga.setHelp( FileManager.readInputStream( SlimDimontPredictorWeb.class.getClassLoader().getResourceAsStream( "projects/slim/helpPredictor.txt" ) ).toString() );//TODO
		
		if(!ga.parse( args, false )){
			System.exit( 1 );
		}	
		
		
		String value = (String)params.getParameterForName( "Value tag" ).getValue();
		String weightingFactor = (String)params.getParameterForName( "Weighting factor" ).getValue();
		double pval = (Double) params.getParameterForName( "p-value" ).getValue();
		
		DataSet data = params.getInputSequences();

		GenDisMixClassifier cl = new GenDisMixClassifier( FileManager.readFile( (String) params.getParameterForName( "SlimDimont" ).getValue() ) );
		
		ThresholdedStrandChIPper model = (ThresholdedStrandChIPper)cl.getDifferentiableSequenceScore( 0 );
		
		Result[][] res = SlimDimontPredictor.run(data,value,weightingFactor,model,pval);
		
		Protocol prot = ga.getProtocol( false );
		
		NumberFormat nf = DecimalFormat.getInstance( Locale.ENGLISH );
		nf.setMaximumFractionDigits( 3 );
		nf.setMinimumFractionDigits( 3 );
		
		for(int i=0;i<res.length;i++){
			prot.append( "\n++++++++++++++++++++++++++++++++++++++++++++++\n\n" );
			prot.appendHeading( "Motif model" );
			prot.append( model.toHtml( nf ) );
			if(res[i].length > 1){
				ga.addResult( res[i][1], false, true);
				ga.addResult( res[i][2], false, true);
				ga.addResult( res[i][3], false, true);
			}
		}
		
		prot.append( "\n++++++++++++++++++++++++++++++++++++++++++++++\n\n" );
		
		prot.appendHeading( "Predicting..." );
		
		for(int i=0;i<res.length;i++){
			ga.addResult( res[i][0], true, false );
			prot.append( "Predicted "+((ListResult)res[i][0]).getRawResult().length+" binding sites." );
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
