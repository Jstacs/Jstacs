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
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.tools.ui.galaxy.DataColumnParameter;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Pair;
import de.jstacs.utils.graphics.GraphicsAdaptorFactory.OutputFormat;



public class DependencyLogoWebParameterSet extends ParameterSet {

	private SimpleParameterSet advanced;
	
	public DependencyLogoWebParameterSet() throws Exception {
		super();
		FileParameter fp = new FileParameter( "Input file", "The tabular input file for plotting the dependency logo", "tabular", true );
		this.parameters.add( fp );
		//this.parameters.add( new SimpleParameter( DataType.INT, "Sequence column", "The index of the column (counting from 1) that contains the sequence data", true, new NumberValidator<Integer>( 1, Integer.MAX_VALUE ), 1 ) );
		this.parameters.add( new DataColumnParameter( fp.getName(), "Sequence column", "The index of the column that contains the sequence data", true ) );
		//this.parameters.add( new SimpleParameter( DataType.INT, "Weight column", "The index of the column (counting from 1) that contains the sequence weights", true, new NumberValidator<Integer>( 1, Integer.MAX_VALUE ), 1 ) );
		this.parameters.add( new DataColumnParameter( fp.getName(), "Weight column", "The index of the column that contains the sequence weights, optional", false ) );
		this.parameters.add( new SelectionParameter( DataType.BOOLEAN, new String[]{"Descending","Ascending"}, new Boolean[]{false,true}, "Order of values", "Sort sequences according to weights ascendingly (smallest value on top) or descendingly (largest values on top). Only has an effect if &quot;Weight column&quot; is set.", false ) );
		
		this.parameters.add( new SelectionParameter( DataType.PARAMETERSET,  
				new String[]{OutputFormat.PDF.name(),OutputFormat.SVG.name(),OutputFormat.PNG.name(),OutputFormat.JPEG.name()},
				new Object[]{new SimpleParameterSet( new SimpleParameter( DataType.INT, "Width", "The width of the output graphic.", true, new NumberValidator<Integer>( 200, Integer.MAX_VALUE ), 2500 ) ),
				             new SimpleParameterSet( new SimpleParameter( DataType.INT, "Width", "The width of the output graphic.", true, new NumberValidator<Integer>( 200, Integer.MAX_VALUE ), 2500 ) ),
				             new SimpleParameterSet( new SimpleParameter( DataType.INT, "Width", "The width of the output graphic (pixels).", true, new NumberValidator<Integer>( 20, Integer.MAX_VALUE ), 500 ) ),
				             new SimpleParameterSet( new SimpleParameter( DataType.INT, "Width", "The width of the output graphic (pixels).", true, new NumberValidator<Integer>( 20, Integer.MAX_VALUE ), 500 ) )},
				"Output format", "The output format of the dependency logo",true
				) );
		
		ParameterList advancedList = new ParameterList();
		
		ExpandableParameterSet exp = new ExpandableParameterSet( getBlockParameters( 250, 300, null ),
				"Numbers of sequences in part", "The number of sequences in the block for part" );
		//exp.addParameterToSet();
		
		advancedList.add( new ParameterSetContainer( "Blocks of sequences","Specify the blocks of sequences that are drawn with their own sequence logo. The remainder of sequences (if enough sequences remaining) is drawn as an additional block.",exp ) );
		
		advancedList.add( new SimpleParameter( DataType.INT, "Height of last block","The height of the last (or only) block of sequences.",true,new NumberValidator<Integer>( 1, Integer.MAX_VALUE ), 750)); 
		
		advancedList.add( new SimpleParameter( DataType.INT, "Height of sequence logo", "The height of the sequence logo(s) of each block", true,new NumberValidator<Integer>( 1, Integer.MAX_VALUE ),250 ) );
		advancedList.add( new SimpleParameter( DataType.INT, "Number of dependencies", "The number of dependencies considered for each position to determine highly dependend positions", true,new NumberValidator<Integer>( 1, Integer.MAX_VALUE ),3 ) );
		
		advanced = new SimpleParameterSet( advancedList.toArray( new Parameter[0] ) );
		
		this.parameters.add( new SelectionParameter( DataType.PARAMETERSET, new String[]{"Default values","Advanced parameters"}, new ParameterSet[]{new SimpleParameterSet( ),advanced}, "Parameters", 
				"If you use default parameters, three blocks of sequences will be drawn, where the first contains 250 sequences, the second contains 1250 sequences, " +
				"and the third contains the remainder of sequences. "+
				"For determining highly dependend positions, the three largest dependencies of each position will be considered.", true ) );
		
	}
	
	private SimpleParameterSet getBlockParameters(int num, int height, String blocknum) throws Exception {
		return new SimpleParameterSet( 
				new SimpleParameter( DataType.INT, "Number of sequences", "The number of sequences in "+(blocknum == null ? "this block" : "block "+blocknum), true,new NumberValidator<Integer>( 1, Integer.MAX_VALUE ),num ), 
				new SimpleParameter( DataType.INT, "Height","The height of "+(blocknum == null ? "this block" : "block "+blocknum)+".",true,new NumberValidator<Integer>( 1, Integer.MAX_VALUE ), height));
	}
	
	public int[] getNumbersOfSequencesForBlocks(){
		if(((SelectionParameter)this.parameters.get( 5 )).getSelected()==0){
			return new int[]{250,1250,0};
		}else{
			ExpandableParameterSet eps = (ExpandableParameterSet)((ParameterSetContainer)advanced.getParameterAt( 0 )).getValue();
			int[] blocks = new int[eps.getNumberOfParameters()+1];
			for(int i=0;i<blocks.length-1;i++){
				blocks[i] = (Integer) ((ParameterSetContainer)eps.getParameterAt( i )).getValue().getParameterAt( 0 ).getValue();
			}
			return blocks;
		}
	}
	
	public int[] getHeightsOfBlocks(){
		if(((SelectionParameter)this.parameters.get( 5 )).getSelected()==0){
			OutputFormat of = getOutputFormat();
			if(of == OutputFormat.PNG || of == OutputFormat.JPEG){
				return new int[]{60,75,150};
			}else{
				return new int[]{300,375,750};
			}
		}else{
			ExpandableParameterSet eps = (ExpandableParameterSet)((ParameterSetContainer)advanced.getParameterAt( 0 )).getValue();
			int[] blocks = new int[eps.getNumberOfParameters()+1];
			for(int i=0;i<blocks.length-1;i++){
				blocks[i] = (Integer) ((ParameterSetContainer)eps.getParameterAt( i )).getValue().getParameterAt( 1 ).getValue();
			}
			blocks[blocks.length-1] = (Integer) advanced.getParameterAt( 1 ).getValue();
			return blocks;
		}
	}
	
	public int getHeightOfSequenceLogo(){
		if(((SelectionParameter)this.parameters.get( 5 )).getSelected()==0){
			return (int)Math.round( getWidth()/(double)10 );
		}else{
			return (Integer)advanced.getParameterAt( 2 ).getValue();
		}
	}
	
	public int getNumberOfDependencies(){
		if(((SelectionParameter)this.parameters.get( 5 )).getSelected()==0){
			return 3;
		}else{
			return (Integer)advanced.getParameterAt( 3 ).getValue();
		}
	}
	
	public int getWidth(){
		SimpleParameterSet val = (SimpleParameterSet)parameters.get( 4 ).getValue();
		return (Integer)val.getParameterAt( 0 ).getValue();
	}
	
	public OutputFormat getOutputFormat(){
		int selected = ((SelectionParameter)parameters.get( 4 )).getSelected();
		String key = ((SelectionParameter)parameters.get( 4 )).getParametersInCollection().getParameterAt( selected ).getName();
		return OutputFormat.valueOf( key );
	}
	
	public Pair<DataSet,double[]> getData() throws Exception {
		int seqCol = (Integer)parameters.get( 1 ).getValue()-1;
		DoubleList w = new DoubleList();;
		int valCol = -1;
		boolean ascending = false;
		if(parameters.get( 2 ).isSet()){
			valCol = (Integer)parameters.get( 2 ).getValue()-1;
			ascending = (Boolean) parameters.get( 3 ).getValue();
		}
		BufferedReader read = new BufferedReader( new FileReader( ( (FileParameter)parameters.get( 0 ) ).getFileContents().getFilename() ) );
		LinkedList<Sequence> seqs = new LinkedList<Sequence>();
		
		String str = null;
		while( (str = read.readLine()) != null ){
			if(!str.trim().startsWith( "#" )){
				String[] parts = str.split( "\t" );

				seqs.add( Sequence.create( DNAAlphabetContainer.SINGLETON, parts[seqCol], "" ) );
				if(valCol > -1){
					w.add( (ascending?-1.0:1.0)*Double.parseDouble( parts[valCol] ) );
				}else{
					w.add( 1.0 );
				}
			}
		}
		read.close();
		
		Sequence[] seqs2 = seqs.toArray(new Sequence[0]);
		double[] w2 = w.toArray();
		
		return new Pair<DataSet, double[]>( new DataSet("",seqs2), w2 );
	}

}
