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

package projects.talen;

import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;

import projects.talen.FastTALENScanner.FastTALENScannerParameterSet;
import projects.talen.FastTALENScanner.Output;
import projects.tals.TALgetter13DiffSM;
import projects.tals.TALgetterDiffSM;
import de.jstacs.DataType;
import de.jstacs.data.GFFParser.GFFList;
import de.jstacs.io.FileManager;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.ParameterSetTagger;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.ListResult;
import de.jstacs.utils.Pair;


public class ScanForTALENCLI {
	
	
	/**
	 * @param args
	 */
	public static void main( String[] args ) throws Exception {
		
		FastTALENScannerParameterSet params = new FastTALENScannerParameterSet();
		
		SelectionParameter mod = new SelectionParameter( DataType.STRING, new String[]{"TALgetter", "TALgetter13"}, new String[]{"TALgetter", "TALgetter13"}, "Model type", 
				"TALgetter is the default model that uses individual binding specificities for each RVD. " +
				"TALgetter13 uses binding specificities that only depend on amino acid 13, i.e., the second amino acid of the repat." +
				"While TALgetter is recommended in most cases, the use of TALgetter13 may be beneficial if you search for target sites of TAL effector with many rare RVDs, for instance YG, HH, or S*.", true );
		mod.setDefault( "TALgetter" );
		
		SelectionParameter filter = (SelectionParameter)params.getParameterForName( "Filter" );
		
		((SelectionParameter)params.getParameterForName( "Architecture" )).setValue( "Custom" );
		
		String filterString = Arrays.toString( FastTALENScanner.filters );
		filterString = filterString.substring( 1, filterString.lastIndexOf( "," ) );
		
		SimpleParameterSet tempPars = new SimpleParameterSet(
				new SimpleParameter( DataType.STRING, "Input sequences", "The sequences to scan for TALEN targets, FastA", true ),
				params.getParameterForName( "Annotation file" ),
				params.getParameterForName( "First RVD sequence" ),
				params.getParameterForName( "Second RVD sequence" ),
				params.getParameterForName( "N-Terminal first" ),
				params.getParameterForName( "N-Terminal second" ),
				params.getParameterForName( "Hetero-dimers only" ),
				((SimpleParameterSet)params.getParameterForName( "Architecture" ).getValue()).getParameterForName( "Minimum distance" ),
				((SimpleParameterSet)params.getParameterForName( "Architecture" ).getValue()).getParameterForName( "Maximum distance" ),
				mod,
				new SimpleParameter( DataType.STRING, "RVD specificities", "File defining additional or overriding existing RVD specificities", false ),
				new SimpleParameter( DataType.DOUBLE, "Filter", filter.getComment()+". Typical values are "+filterString, true, new NumberValidator<Double>( 0.35, 1.0 ), 0.4 ),
				params.getParameterForName( "Maximum number of targets" ),
				new SimpleParameter( DataType.STRING, "Additional output", "Path to a GFF3/GFF2 file to which predictions are written in addition to the default output, extension defines format (.gff3/.gff)", false ),
				new SimpleParameter( DataType.INT, "Number of threads", "Number of threads used by TALoffer. More than 3 threads typically do not lead to an additional speed-up.", true, new NumberValidator<Integer>(1,8), 3 )
				);
		
		ParameterSetTagger pst = new ParameterSetTagger( 
				new String[]{"input","annotation", "rvdl" ,"rvdr","nterml", "ntermr","heterodimers","min", "max", "model","addrvds","filter","top","out","numThreads"}, 
				tempPars);
		
		pst.fillParameters( "=", args );
		
		if(tempPars.getParameterForName( "RVD specificities" ).isSet()){
			params.setRVDSpecificities((String)tempPars.getParameterForName( "RVD specificities" ).getValue());
		}
		
		if(tempPars.getParameterForName( "Additional output" ).isSet() ){
			String file = (String)tempPars.getParameterForName( "Additional output" ).getValue();
			if(file.toLowerCase().endsWith( "gff" )){
				params.getParameterForName( Output.class.getSimpleName() ).setValue( Output.GFF2 );
			}else{
				params.getParameterForName( Output.class.getSimpleName() ).setValue( Output.GFF3 );
			}
		}else{
			params.getParameterForName( Output.class.getSimpleName() ).setValue( Output.NONE );
		}
		
		Double value = (Double)tempPars.getParameterForName( "Filter" ).getValue();
		params.getParameterForName( "Filter" ).setValue( "Custom" );
		((SimpleParameterSet)params.getParameterForName( "Filter" ).getValue()).getParameterAt( 0 ).setValue( value );
		
		
		System.err.println( "parameters:" );
		System.err.println( pst );
		System.err.println("_________________________________");
		if( !tempPars.hasDefaultOrIsSet() ) {
			System.err.println( "Some of the required parameters are not specified." );
			System.exit( 1 );
		}
		params.setInputPath( (String)tempPars.getParameterAt( 0 ).getValue() );
		
		TALgetterDiffSM model = null;
		
		if(mod.getValue().equals( "TALgetter" )){
			model = (TALgetterDiffSM)XMLParser.extractObjectForTags( FileManager.readInputStream( ScanForTALENCLI.class.getClassLoader().getResourceAsStream( "projects/tals/talfinder_obg2_hyp_bg.xml" ) ), "model" );		
		}else{
			model = (TALgetter13DiffSM)XMLParser.extractObjectForTags( FileManager.readInputStream( ScanForTALENCLI.class.getClassLoader().getResourceAsStream( "projects/tals/talfinder_obg2_hyp_bg_map.xml" ) ), "model" );
		}
		
		Pair<String[], double[][][]> specs = params.getSpecificities();
		if(specs != null){
			model.addAndSet( specs.getFirstElement(), specs.getSecondElement()[0], specs.getSecondElement()[1][0] );
		}
		
		
		((TALgetterDiffSM)model).fix();
		
		//System.out.println(model);
		
		BufferedReader ds = params.getInputSequences();
		
		System.err.println( "Searching for targets of TALEN" );
		System.err.println( "\t"+params.getLeftTALSequence()+" & "+params.getRightTALSequence()+" with distance between "+params.getMinimumDistance()+" and "+params.getMaximumDistance()+".\n" );
		System.err.println( "Reporting best "+params.getN()+" target sites.\n\n" );
		
		GFFList gff = null;
		try{
			gff = params.getAnnotation();
		}catch(Exception e){
			gff = null;
			System.err.println("Could not load Annotation.");
		}
		
		ListResult[] rls = (new FastTALENScanner()).scan( (TALgetterDiffSM)model, params, ds, gff, (Integer)tempPars.getParameterForName( "Number of threads" ).getValue() );
		
		
		System.out.println( rls[0] );

		System.out.flush();
		
		if(tempPars.getParameterForName( "Additional output" ).isSet() ){
			String file = (String)tempPars.getParameterForName( "Additional output" ).getValue();
			File file2 = new File(file);
			if(file2.exists()){
				System.err.println("Additional output file already exists. Please choose different file.");
			}else{
				PrintWriter wr = new PrintWriter( file );
				if(file.toLowerCase().endsWith( "gff" )){
					wr.println("##gff-version 2");
				}else{
					wr.println("##gff-version 3");
				}
				wr.println(rls[1]);
				wr.close();
			}
		}
		//System.err.println("TIME: "+InfixMatchFinder.time);
		System.err.println( "Finished..." );
		
	}
	
	
}
