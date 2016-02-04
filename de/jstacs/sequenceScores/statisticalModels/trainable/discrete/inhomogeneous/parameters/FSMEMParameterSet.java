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

package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters;

import de.jstacs.DataType;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.ConstraintManager.Decomposition;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.FSMEManager;

/**
 * The ParameterSet for a FSMEManager.
 * 
 * @author Jens Keilwagen
 * 
 * @see FSMEManager
 */
public class FSMEMParameterSet extends MEManagerParameterSet
{
	private static final String[] structureNames = new String[]{ "m2sx", "m3sx", "m1s0", "m2s0", "m3s0", "m4s0",
			"m5s0", "alternative" };

	private static final String[] structureComments = new String[]{ "all dependencies between 2 nodes will be used",
			"all dependencies between 3 nodes", "no dependencies will be used",
			"only nearest neighbour dependencies will be used",
			"nearest neighbour and 2nd nearest neighbour dependencies will be used",
			"1st, 2nd and 3rd nearest neighbour dependencies will be used",
			"1st, 2nd, 3rd and 4th nearest neighbour dependencies will be used",
			"here you can input your structure by your own (nearly everything is possible)" };

	private static final Object[] structure = new Object[]{ "m2sx", "m3sx", "m1s0", "m2s0", "m3s0", "m4s0", "m5s0",
			new ConstraintParameterSet() };

	/**
	 * The simple constructor.
	 * 
	 * @throws Exception if some parameters templates could not be created properly
	 */
	public FSMEMParameterSet() throws Exception
	{
		super( FSMEManager.class );
	}

	/**
	 * The constructor for the {@link de.jstacs.Storable} interface.
	 * 
	 * @param s
	 *            the StringBuffer
	 * 
	 * @throws NonParsableException
	 *             if the StringBuffer is not parsable
	 */
	public FSMEMParameterSet( StringBuffer s ) throws NonParsableException
	{
		super( s );
	}

	/**
	 * The fast constructor.
	 * 
	 * @param alphabet
	 *            the alphabet
	 * @param length
	 *            the length of the modeled sequences
	 * @param ess
	 *            the ess
	 * @param description
	 *            the description
	 * @param decomposition
	 *            the kind of decomposition
	 * @param reduce
	 *            whether the constraints should be reduced
	 * @param algorithm
	 *            the choice of algorithm
	 * @param epsilon
	 *            the threshold for stopping the numerical algorithms
	 * @param constraints
	 *            the constraints to be used<br>
	 *            There are two ways entering the constraints.
	 *            <ol>
	 *            <li> Enter <b>single constraints</b> in the following way: node_1, node_2, ..., node_k.
	 *            <li> Enter <b>groups of constraints</b> in the following way: &quot;m&lt;number of
	 *            nodes&gt;s&lt;number of skipped nodes&gt;&quot;. E.g. m2s1 is the set {0,2; 1,3; 2,4; ...}.<br>
	 *            The &quot;x&quot; can be used to determine that each skip should be used (e.g. m2sx).
	 *            </ol>
	 *            To separate the constraints use &quot;;&quot;.<br>
	 *            <br>
	 *            Additionally for each node a constraint is assumed (i.e. {0; 1; 2; ...}).
	 * 
	 * @throws Exception
	 *             if something went wrong
	 *             
	 * @see Decomposition#DECOMPOSE_NOTHING
	 * @see Decomposition#DECOMPOSE_UNCONNECTED
	 * @see Decomposition#DECOMPOSE_LESS_CONNECTED
	 */
	public FSMEMParameterSet( AlphabetContainer alphabet, int length, double ess, String description,
			Decomposition decomposition, boolean reduce, byte algorithm, double epsilon, String constraints ) throws Exception
	{
		super( FSMEManager.class, alphabet, length, ess, description, decomposition, reduce, algorithm, epsilon );
		int index = getIndex( structureNames, structure, constraints, true );
		parameters.get( 6 ).setValue( structureNames[index] );
		if( index == structureNames.length - 1 )
		{
			((ParameterSet) parameters.get( 6 ).getValue()).getParameterAt( 0 ).setValue( constraints );
		}
	}

	protected void addParameters() throws Exception
	{
		super.addParameters();
		parameters
				.add( new SelectionParameter(
						DataType.STRING,
						structureNames,
						structure,
						structureComments,
						"structure",
						"The structure is very important for defining a probability distribustion. It states which variables dependent on each other.",
						true ) );
		parameters.get( 6 ).setDefault( structureNames[0] );
	}

	public String getInstanceComment()
	{
		return "holds the parameters for a CMEManager";
	}
}