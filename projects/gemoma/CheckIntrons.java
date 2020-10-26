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

package projects.gemoma;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import de.jstacs.DataType;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;

/**
 * The class creates statistics for introns.
 *  
 * @author Jens Keilwagen
 */
public class CheckIntrons extends GeMoMaModule {

	@Override
	public ToolParameterSet getToolParameters() {
		try{
			return new ToolParameterSet( getShortName(), 
					new FileParameter( "target genome", "The target genome file (FASTA), i.e., the target sequences in the blast run. Should be in IUPAC code", "fasta", true ),					
					new ParameterSetContainer( "introns", "", new ExpandableParameterSet( new SimpleParameterSet(	
							new FileParameter( "introns", "Introns (GFF), which might be obtained from RNA-seq", "gff", true )
						), "introns", "", 1 ) ),
					
					new SimpleParameter( DataType.BOOLEAN, "verbose", "A flag which allows to output a wealth of additional information per transcript", true, false )
			);
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, String temp) throws Exception {
		boolean verbose = (Boolean) parameters.getParameterForName("verbose").getValue();
		String targetGenome = (String) parameters.getParameterForName("target genome").getValue(); 
		int reads = 1;
		ExpandableParameterSet introns = (ExpandableParameterSet)((ParameterSetContainer)parameters.getParameterAt(1)).getValue();
		
		HashMap<String,String> seqs = Tools.getFasta(targetGenome,20);
		
		ArrayList<String> fName = new ArrayList<String>();
		for( int i = 0; i < introns.getNumberOfParameters(); i++ ) {
			Parameter y = ((ParameterSet)introns.getParameterAt(i).getValue()).getParameterAt(0);
			if( y.isSet() ) {
				fName.add(y.getValue().toString());
			}
		}
		if( fName.size()>0 ) {
			HashMap<String,int[]> diNucl = new HashMap<String,int[]>();
			HashMap<String, int[][][]>[] res = GeMoMa.readIntrons( reads, protocol, verbose, seqs, diNucl, fName.toArray(new String[fName.size()]) );
			String[] keys = diNucl.keySet().toArray(new String[0]);
			Arrays.sort( keys );
			for( int i = 0; i < keys.length; i++ ) {
				int[] v = diNucl.get( keys[i] );
				protocol.append( keys[i] + "\t" + v[0] + "\t" + v[1] + "\n" );
			}
		}
		return null;
	}

	@Override
	public String getToolName() {
		return "CheckIntrons";
	}

	@Override
	public String getShortName() {
		return getToolName();
	}

	@Override
	public String getDescription() {
		return "creates statistics for introns";
	}

	@Override
	public String getHelpText() {
		return "The tool checks the distribution of introns on the strands and the dinucleotide distribution at splice sites." + MORE;
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}
	
	@Override
	public ToolResult[] getTestCases( String path ) {
		// TODO missing test cases
		return null;
	}
}
