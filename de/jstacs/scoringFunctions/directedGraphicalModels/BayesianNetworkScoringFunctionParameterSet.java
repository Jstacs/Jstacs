/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.scoringFunctions.directedGraphicalModels;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.AlphabetContainer.AlphabetContainerType;
import de.jstacs.io.ParameterSetParser.NotInstantiableException;
import de.jstacs.parameters.InstanceParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SequenceScoringParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.measures.Measure;
import de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.measures.Measure.MeasureParameterSet;
import de.jstacs.utils.SubclassFinder;

/**
 * Class for the parameters of a {@link BayesianNetworkScoringFunction}. This
 * class fulfills the requirements of a {@link SequenceScoringParameterSet} and
 * can be used to create a new {@link BayesianNetworkScoringFunction}.
 * 
 * @author Jan Grau
 */
public class BayesianNetworkScoringFunctionParameterSet extends
		SequenceScoringParameterSet {

	/**
	 * Creates a new {@link BayesianNetworkScoringFunctionParameterSet} with
	 * pre-defined parameter values.
	 * 
	 * @param alphabet
	 *            the alphabet of the scoring function boxed in an
	 *            {@link AlphabetContainer}, e.g
	 *            <code>new AlphabetContainer(new DNAAlphabet())</code>
	 * @param length
	 *            the length of the scoring function, i.e. the length of the
	 *            sequences this scoring function can handle
	 * @param ess
	 *            the equivalent sample size
	 * @param plugInParameters
	 *            indicates if plug-in parameters, i.e. generative (MAP)
	 *            parameters, shall be used upon initialization
	 * @param structureMeasure
	 *            the {@link Measure} used for the structure, e.g.
	 *            {@link de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.measures.InhomogeneousMarkov}
	 * 
	 * @throws Exception
	 *             if the alphabet or the length are not in the expected range
	 *             of values
	 */
	public BayesianNetworkScoringFunctionParameterSet(
			AlphabetContainer alphabet, int length, double ess,
			boolean plugInParameters, Measure structureMeasure)
			throws Exception {
		super(BayesianNetworkScoringFunction.class, alphabet, length, false);
		addParameters();
		parameters.get(0).setValue(ess);
		parameters.get(1).setValue(plugInParameters);
		InstanceParameterSet struct = structureMeasure.getCurrentParameterSet();
		parameters.get(2).setValue(struct);
		((ParameterSetContainer) ((SelectionParameter) parameters.get(2))
				.getParametersInCollection()
				.getParameterAt(
						((SelectionParameter) parameters.get(2)).getSelected()))
				.setValue(struct);
	}

	/**
	 * Creates a new {@link BayesianNetworkScoringFunctionParameterSet} with
	 * empty parameter values.
	 * @throws Exception 
	 */
	public BayesianNetworkScoringFunctionParameterSet() throws Exception {
		super(BayesianNetworkScoringFunction.class,
				AlphabetContainerType.DISCRETE, false);
		addParameters();
		
	}

	private void addParameters() throws Exception {
		parameters.add(new SimpleParameter(DataType.DOUBLE, "ESS",
				"The equivalent sample size", true));
		parameters.add(new SimpleParameter(DataType.BOOLEAN,
				"Plug-in parameters", "Use plug-in parameters", true));
		parameters.add(SubclassFinder.getCollection(MeasureParameterSet.class,
				Measure.class.getPackage().getName(), "Structure measure",
				"Choose a measure to determine the structure.", true));
	}
	
	/**
	 * Creates a new {@link BayesianNetworkScoringFunctionParameterSet} from its
	 * XML representation as defined by the {@link de.jstacs.Storable}
	 * interface.
	 * 
	 * @param representation
	 *            the XML code as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             is thrown if the XML representation could not be parsed
	 */
	public BayesianNetworkScoringFunctionParameterSet(
			StringBuffer representation) throws NonParsableException {
		super(representation);
	}

	/**
	 * Returns the equivalent samples size (ess) defined in this set of
	 * parameters.
	 * 
	 * @return the ess
	 */
	public double getEss() {
		return (Double) parameters.get(0).getValue();
	}

	/**
	 * Returns true if plug-in parameters shall be used when creating a
	 * {@link BayesianNetworkScoringFunction} from this set of parameters.
	 * 
	 * @return if plug-in parameters shall be used
	 */
	public boolean getPlugInParameters() {
		return (Boolean) parameters.get(1).getValue();
	}

	/**
	 * Returns the structure {@link Measure} defined by this set of parameters.
	 * 
	 * @return the structure {@link Measure}
	 * 
	 * @throws NotInstantiableException
	 *             if the {@link Measure} could not be created from its own
	 *             {@link InstanceParameterSet}
	 */
	public Measure getMeasure() throws NotInstantiableException {
		return (Measure) ((InstanceParameterSet) parameters.get(2).getValue())
				.getInstance();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
	 */
	@Override
	public String getInstanceComment() {
		return "Scoring function for all special cases of moral Bayesian networks";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "Bayesian network scoring function";
	}

}
