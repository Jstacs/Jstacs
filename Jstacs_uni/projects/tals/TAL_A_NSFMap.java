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

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.DiscreteAlphabetMapping;
import de.jstacs.data.alphabets.DoubleSymbolException;
import de.jstacs.io.NonParsableException;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.utils.ToolBox;

/**
 * Sub-class of {@link TALgetterRVDDependentComponent}, which assigns the same specificities to all RVDs with the same 13th amino acid.
 * Used in {@link TALgetter13DiffSM}.
 * 
 * @author Jan Grau
 *
 */
public class TAL_A_NSFMap extends TALgetterRVDDependentComponent {

	private int[] map;
	
	/**
	 * Creates a new {@link TAL_A_NSFMap} object for given alphabet and RVD alphabet, prior importances and prior 
	 * binding preferences.
	 * 
	 * @param alphabets the alphabet, typically {@link DNAAlphabet}
	 * @param alphabetsRVD the RVD alphabet
	 * @param length the expected length of target sites
	 * @param ess the equivalent sample size of the prior
	 * @param priorImp the prior importances, in same order as RVDs in <code>alphabetsRVD</code>
	 * @param priorPrefs the prior binding preferences, in same order as RVDs in <code>alphabetsRVD</code>
	 * @throws Exception is something went wrong
	 */
	public TAL_A_NSFMap( AlphabetContainer alphabets, AlphabetContainer alphabetsRVD, int length, double ess, double[] priorImp,
							double[][] priorPrefs ) throws Exception {
		super( alphabets, alphabetsRVD, length, ess, priorImp, priorPrefs );
	}
	
	
	
	@Override
	public TAL_A_NSFMap clone() throws CloneNotSupportedException {
		TAL_A_NSFMap clone = (TAL_A_NSFMap)super.clone();
		if(map != null){
			clone.map = map.clone();
		}
		return clone;
	}



	private static int[] buildMap(AlphabetContainer con){
		HashSet<String> set = new HashSet<String>();
		for(int i=0;i<con.getAlphabetLengthAt( 0 );i++){
			String sym = con.getSymbol( 0, i ).toUpperCase();
			String last = sym.substring( 1 );
			set.add( last );
		}
		try{
			DiscreteAlphabet da = new DiscreteAlphabet( true, set.toArray( new String[0] ) );
			
			int[] map = new int[(int)con.getAlphabetLengthAt( 0 )];
			
			for(int i=0;i<con.getAlphabetLengthAt( 0 );i++){
				String sym = con.getSymbol( 0, i ).toUpperCase();
				String last = sym.substring( 1 );
				map[i] = da.getCode( last );
			}
			return map;
		}catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}

	@Override
	protected int getNumberOfSymbols( AlphabetContainer con ) {
		if(map == null){
			map = buildMap(con);
		}
		int max = -1;
		for(int i=0;i<map.length;i++){
			if(map[i] > max){
				max = map[i];
			}
		}
		return max+1;
	}

	@Override
	protected int getMappedIndex( AlphabetContainer con, int original ) {
		if(map == null){
			map = buildMap(con);
		}
		return map[original];
	}
	
	@Override
	public AlphabetContainer addAndSet(String[] rvds, double[][] specs) throws WrongAlphabetException, IllegalArgumentException, DoubleSymbolException{
		DiscreteAlphabet alph = (DiscreteAlphabet)alphabetsRVD.getAlphabetAt( 0 );
		int[] map = new int[rvds.length];
		Arrays.fill( map, -1 );
		int n = 0;
		for(int i=0;i<rvds.length;i++){
			if(alph.isSymbol( rvds[i] )){
				map[i] = alph.getCode( rvds[i] );
			}else{
				n++;
			}
		}
		if(n > 0){
			String[] newAlph = new String[(int)alph.length()+n];
			int i=0;
			for(;i<alph.length();i++){
				newAlph[i] = alph.getSymbolAt( i );
			}
			for(int j=0;j<map.length;j++){
				if(map[j]==-1){
					newAlph[i] = rvds[j];
					i++;
				}
			}
			alph = new DiscreteAlphabet( true, newAlph );
			alphabetsRVD = new AlphabetContainer( alph );
			this.map = buildMap( alphabetsRVD );
			
			n = getNumberOfSymbols( alphabetsRVD );
			HomogeneousMMDiffSM[] temp_c = new HomogeneousMMDiffSM[n];
			System.arraycopy( hmm_c, 0, temp_c, 0, hmm_c.length );
			for(i=hmm_c.length;i<n;i++){
				temp_c[i] = new HomogeneousMMDiffSM(alphabets,0,ess/temp_c.length,priorLength);
				temp_c[i].initializeFunctionRandomly( false );
			}
			hmm_c = temp_c;
		}
		
		for(int i=0;i<rvds.length;i++){
			int idx =getMappedIndex( alphabetsRVD, alph.getCode( rvds[i] ) );
			double[] temp = specs[i].clone();
			for(int j=0;j<temp.length;j++){
				temp[j] = Math.log( temp[j] );
			}
			hmm_c[idx].setParameters( temp, 0 );
		}
		
		return alphabetsRVD;		
	}

	/**
	 * Creates a new {@link TAL_A_NSFMap} from its XML description
	 * @param xml the XML description
	 * @throws NonParsableException if the XML could not be parsed
	 */
	public TAL_A_NSFMap( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	
	@Override
	public String toString(NumberFormat nf){
		StringBuffer sb = new StringBuffer();
		for(int i=0;i<alphabetsRVD.getAlphabetLengthAt( 0 );i++){
			sb.append( alphabetsRVD.getSymbol( 0, i )+"\t" );
			sb.append( hmm_c[getMappedIndex( alphabetsRVD, i )]+"\n" );
		}
		return sb.toString();
	}
	
}
