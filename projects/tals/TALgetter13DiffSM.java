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

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.io.NonParsableException;

/**
 * Class for the TALgetter13 variant of TALgetter.
 * 
 * @author Jan Grau
 *
 */
public class TALgetter13DiffSM extends TALgetterDiffSM {

	/**
	 * Creates a new TALgetter13 model.
	 *  
	 * @param alphabets the alphabet of target sites, typically {@link DNAAlphabetContainer}.
	 * @param alphabetsRVD the alphabet of RVDs
	 * @param midLength the expected length of target sites
	 * @param ess the equivalent sample size
	 * @param order_talU the order of the RVD-independent component
	 * @param priorFP the prior probabilities for position 0
	 * @param priorImp the prior importances, in same order as RVDs in <code>alphabetsRVD</code>
	 * @param priorPrefs the prior binding preferences, in same order as RVDs in <code>alphabetsRVD</code>
	 * @throws Exception if something went wrong
	 */
	public TALgetter13DiffSM( AlphabetContainer alphabets, AlphabetContainer alphabetsRVD, double midLength, double Ess,
								int order_talU,double[] priorFP, double[] priorImp, double[][] priorPrefs ) throws Exception {
		super( alphabets, alphabetsRVD, midLength, Ess, order_talU, priorFP, priorImp, priorPrefs );
	}

	/**
	 * Creates a new {@link TALgetter13DiffSM} from its XML description.
	 * @param xml the XML description
	 * @throws NonParsableException if the description could not be parsed.
	 */
	public TALgetter13DiffSM( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	protected TALgetterRVDDependentComponent getTalANsf( AlphabetContainer alphabets, AlphabetContainer alphabetsRVD, int midLength, double Ess, double part,
			double[] priorImp, double[][] priorPrefs ) throws Exception {
		return new TAL_A_NSFMap(alphabets, alphabetsRVD, (int)midLength, Ess*part, priorImp,priorPrefs);
	}

	
	
}
