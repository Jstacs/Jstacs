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

import java.awt.image.BufferedImage;
import java.io.BufferedReader;

import de.jstacs.DataType;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.FileManager;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ImageResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor.Protocol;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SeqLogoPlotter;
import projects.talen.FastTALENScanner.FastTALENScannerParameterSet;
import projects.talen.FastTALENScanner.Output;
import projects.talen.GFFParser.GFFList;
import projects.tals.ScanForTBSWeb;
import projects.tals.TALgetter13DiffSM;
import projects.tals.TALgetterDiffSM;


public class ScanForTALENWeb {
	
	
	/**
	 * @param args
	 */
	public static void main( String[] args ) throws Exception {
		
		FastTALENScannerParameterSet params = new FastTALENScannerParameterSet();
		
		SelectionParameter mod2 = new SelectionParameter( DataType.STRING, new String[]{"TALgetter", "TALgetter13"}, new String[]{"TALgetter", "TALgetter13"}, "Model type", 
				"TALgetter is the default model that uses individual binding specificities for each RVD. " +
				"TALgetter13 uses binding specificities that only depend on amino acid 13, i.e., the second amino acid of the repat." +
				"While TALgetter is recommended in most cases, the use of TALgetter13 may be beneficial if you search for target sites of TALENs with many rare RVDs, for instance YG, HH, or S*.", true );
		mod2.setDefault( "TALgetter" );

		
		SelectionParameter train = new SelectionParameter( DataType.PARAMETERSET, new String[]{"Use standard model", "Use previously trained model"}, new ParameterSet[]{
					                                                                                                                                         new SimpleParameterSet( mod2 ),
					                                                                                                                                         new SimpleParameterSet( new FileParameter( "Model", "Choose a TALgetter model from your history", "xml", true ) )
					}, "Model", "You can either use the standard TALgetter model or re-use a TALgetter model that has already been trained on given pairs of TAL effectors and target sites using the TALgetter application. ", true );
		
		
		SimpleParameterSet tempPars = new SimpleParameterSet(
				params.getParameterForName( "Input sequences" ),
				params.getParameterForName( "Annotation file" ),
				params.getParameterForName( "First RVD sequence" ),
				params.getParameterForName( "N-Terminal first" ),
				params.getParameterForName( "Second RVD sequence" ),
				params.getParameterForName( "N-Terminal second" ),
				params.getParameterForName( "Hetero-dimers only" ),
				params.getParameterForName( "Architecture" ),
				params.getParameterForName( "Filter" ),
				train,
				params.getParameterForName( "RVD specificities" ),
				params.getParameterForName( "Maximum number of targets" ),
				params.getParameterForName( Output.class.getSimpleName() )
				);
		
		
		boolean[] line = new boolean[]{true,false,true,false,true,false,true,true,true,true,false,true,true};
		
		
		GalaxyAdaptor ga = new GalaxyAdaptor( tempPars,null, line,"TALENoffer", "TALENoffer is a tool for predicting off-targets of TAL effector nucleases (TALENs).", "1.0", "java -Xms256M -Xmx4G -jar "+System.getProperty( "user.dir" )+System.getProperty( "file.separator" )+"TALENofferWeb.jar", "jobname", true );
		ga.setHelp( FileManager.readInputStream( ScanForTBSWeb.class.getClassLoader().getResourceAsStream( "projects/talen/helpoffer.txt" ) ).toString() );
		
		if(!ga.parse( args, false )){
			System.exit( 1 );
		}	
		
		
		TALgetterDiffSM model = null;
		
		if(((SelectionParameter)tempPars.getParameterForName( "Model" )).getSelected() == 1){
			FileParameter fp = (FileParameter)((SimpleParameterSet)tempPars.getParameterForName( "Model" ).getValue()).getParameterAt( 0 );
			
			String fn = (String)fp.getValue();
			StringBuffer sb = FileManager.readFile( fn );
			if(mod2.getValue().equals( "TALgetter" )){
				model = new TALgetterDiffSM( sb );
			}else{
				model = new TALgetter13DiffSM( sb );
			}
		}else{

			if(mod2.getValue().equals( "TALgetter" )){
				model = (TALgetterDiffSM)XMLParser.extractObjectForTags( FileManager.readInputStream( ScanForTALENWeb.class.getClassLoader().getResourceAsStream( "projects/tals/talfinder_obg2_hyp_bg.xml" ) ), "model" );		
			}else{
				model = (TALgetterDiffSM)XMLParser.extractObjectForTags( FileManager.readInputStream( ScanForTALENWeb.class.getClassLoader().getResourceAsStream( "projects/tals/talfinder_obg2_hyp_bg_map.xml" ) ), "model" );
			}
		}
		
		Pair<String[], double[][][]> specs = params.getSpecificities();
		if(specs != null){
			model.addAndSet( specs.getFirstElement(), specs.getSecondElement()[0], specs.getSecondElement()[1][0] );
		}
		
		((TALgetterDiffSM)model).fix();
		
		
		
		BufferedReader ds = params.getInputSequences();
		
		Protocol prot = ga.getProtocol( false );
		
		prot.append( "Searching for targets of TALEN<br />" );
		prot.append( "\t"+params.getLeftTALSequence()+" and "+params.getRightTALSequence()+" with distance between "+params.getMinimumDistance()+" and "+params.getMaximumDistance()+".<br />" );
		prot.append( "Reporting best "+params.getN()+" target sites.<br /><br />" );
		
		GFFList gff = null;
		try{
			gff = params.getAnnotation();
		}catch(Exception e){
			gff = null;
			prot.append("Could not load Annotation.<br />");
		}
		
		Pair<double[][],double[]> specNImp = ((TALgetterDiffSM)model).getSpecificitiesAndImportances( Sequence.create( ((TALgetterDiffSM)model).getRVDAlphabet(), params.getLeftTALSequence(), "-" ) );
		
		int height = SeqLogoPlotter.getHeight( 750, specNImp.getFirstElement() );
		
		BufferedImage img = SeqLogoPlotter.plotTALgetterLogoToBufferedImage( height, specNImp.getFirstElement(), specNImp.getSecondElement(), ("0-"+params.getLeftTALSequence()).split( "-" ) );
		
		ga.addResult( new ImageResult( "First monomer target site logo", "Logo plot representing specificities and importances for the first TALEN monomer according to model parameters", img ), false, true );
	
		specNImp = ((TALgetterDiffSM)model).getSpecificitiesAndImportances( Sequence.create( ((TALgetterDiffSM)model).getRVDAlphabet(), params.getRightTALSequence(), "-" ) );
		
		height = SeqLogoPlotter.getHeight( 750, specNImp.getFirstElement() );
		
		img = SeqLogoPlotter.plotTALgetterLogoToBufferedImage( height, specNImp.getFirstElement(), specNImp.getSecondElement(), ("0-"+params.getRightTALSequence()).split( "-" ) );
		
		ga.addResult( new ImageResult( "Second monomer target site logo", "Logo plot representing specificities and importances for the second TALEN monomer according to model parameters", img ), false, true );
	
		
		ListResult[] rls = (new FastTALENScanner()).scan( (TALgetterDiffSM)model, params, ds, gff, 1 );
		
		
		
		ga.addResult( new ListResult( "Description of output columns", "The output can also be downloaded as a tab-separated file.", null, 
				new ResultSet[]{
				                new ResultSet( new Result[]{new CategoricalResult("Column","","ID"),
				                                            new CategoricalResult( "Description", "", "The ID of the input sequence as given in FastA header" )} ),
	                            new ResultSet( new Result[]{new CategoricalResult("Column","","Position 1"),
	                                                        new CategoricalResult( "Description", "", "Position of the first TALEN monomer target site in the given sequence" )} ),
	                            new ResultSet( new Result[]{new CategoricalResult("Column","","Position 2"),
	                        	                            new CategoricalResult( "Description", "", "Position of the second TALEN monomer target site in the given sequence" )} ),
				                new ResultSet( new Result[]{new CategoricalResult("Column","","Distance"),
				                                            new CategoricalResult( "Description", "", "The distance between the two TALEN monomer target sites" )} ),
				                new ResultSet( new Result[]{new CategoricalResult("Column","","Sequence 1"),
				            				                new CategoricalResult( "Description", "", "The sequence of the target site of the first TALEN monomer" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Matches 1"),
				            				                new CategoricalResult( "Description", "", "The matches of the first TALEN monomer. Categories of matches per position: M/m first position match/mismatch; | match; : weak match; x mismatch" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Sequence 2"),
				            				            	new CategoricalResult( "Description", "", "The sequence of the target site of the second TALEN monomer" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Matches 2"),
				            				            	new CategoricalResult( "Description", "", "The matches of the second TALEN monomer. Categories of matches per position: M/m first position match/mismatch; | match; : weak match; x mismatch" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Architecture"),
				            				            	new CategoricalResult( "Description", "", "The architecture of the TALEN dimer responsible for the off-target" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Full site"),
				            				            	new CategoricalResult( "Description", "", "The sequence of both sites and enclosed linker on the forward strand" )} ),				            	
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Score"),
				            				                new CategoricalResult( "Description", "", "Relative score returned by the model" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Features"),
				            				                new CategoricalResult( "Description", "", "Sequence features overlapping the off-target (if annotation file provided)" )} )
				            				                                            
				            				                                            
				            				                                            
				            				                                            
				}
		), false, true );
		
		ga.addResult( rls[0], true, true );
		if(params.getOutput() != Output.NONE){
			ga.addResult( rls[1], true, false, params.getOutput() == Output.GFF3 ? "gff3" : "gff" );
		}
		
		prot.append( "Finished..." );
		
		ga.writeOutput();
		
	}
	
	
}
