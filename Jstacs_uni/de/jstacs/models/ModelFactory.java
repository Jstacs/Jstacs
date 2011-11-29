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

package de.jstacs.models;

import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.io.ParameterSetParser.NotInstantiableException;
import de.jstacs.models.discrete.homogeneous.HomogeneousMM;
import de.jstacs.models.discrete.homogeneous.parameters.HomMMParameterSet;
import de.jstacs.models.discrete.inhomogeneous.BayesianNetworkModel;
import de.jstacs.models.discrete.inhomogeneous.FSDAGModel;
import de.jstacs.models.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.models.discrete.inhomogeneous.StructureLearner.ModelType;
import de.jstacs.models.discrete.inhomogeneous.parameters.BayesianNetworkModelParameterSet;
import de.jstacs.models.discrete.inhomogeneous.parameters.FSDAGMParameterSet;
import de.jstacs.models.mixture.MixtureModel;
import de.jstacs.models.mixture.StrandModel;
import de.jstacs.models.mixture.AbstractMixtureModel.Parameterization;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import de.jstacs.parameters.SequenceScoringParameterSet;


/**
 * This class allows to easily create some frequently used models.
 * It offers only one way of creating each model and set some of the parameters to default values.
 * If you like to set further models please check the constructors of the individual classes. 
 * 
 * @author Jens Keilwagen
 */
public class ModelFactory {

	private static AbstractModel getModelFromParameterSet( SequenceScoringParameterSet ps ) throws NotInstantiableException {
		return (AbstractModel) ps.getInstance();
	}
	
	/**
	 * This method returns a position weight matrix (PWM). A PWM assumes that all positions of a sequence are statistically independent.
	 * 
	 * @param con the {@link AlphabetContainer} of the PWM
	 * @param length the length of the PWM, i.e., the length of the sequences that can be modeled
	 * @param ess the equivalent sample size (ess) of the PWM, if 0 (zero) the model can be trained using maximum likelihood principle, otherwise it can be trained using the maximum a posteriori principle using the BDeu prior
	 * 
	 * @return the PWM
	 * 
	 * @throws Exception if the model can not be created correctly
	 */
	public static FSDAGModel createPWM( AlphabetContainer con, int length, double ess ) throws Exception {
		FSDAGMParameterSet ps = new FSDAGMParameterSet( con, length, ess, null, "" );
		return (FSDAGModel) getModelFromParameterSet( ps );
	}
	
	private static BayesianNetworkModel createBN( AlphabetContainer con, int length, double ess, ModelType type, byte order ) throws Exception {
		BayesianNetworkModelParameterSet ps = new BayesianNetworkModelParameterSet( con, length, ess, null, type, order, LearningType.ML_OR_MAP );
		return (BayesianNetworkModel) getModelFromParameterSet( ps );
	}

	/**
	 * This method returns a inhomogeneous Markov model (IMM) with user-specified order.
	 * 
	 * @param con the {@link AlphabetContainer} of the IMM
	 * @param length the length of the IMM, i.e., the length of the sequences that can be modeled
	 * @param ess the equivalent sample size (ess) of the IMM, if 0 (zero) the model can be trained using maximum likelihood principle, otherwise it can be trained using the maximum a posteriori principle using the BDeu prior
	 * @param order the order of the IMM, i.e., the number of directly preceding random variables (=positions) that might have an influence on the probability of outcome of a random variable (=position) 
	 * 
	 * @return the IMM
	 * 
	 * @throws Exception if the model can not be created correctly
	 */
	public static BayesianNetworkModel createInhomogeneousMarkovModel( AlphabetContainer con, int length, double ess, byte order ) throws Exception {
		return createBN( con, length, ess, ModelType.IMM, order );
	}
	
	/**
	 * This method returns a permuted Markov model (PMM) with user-specified order. Permuted Markov models determine a permutation of the random variables (=position) and than build a inhomogeneous Markov model based on this permutation.
	 * 
	 * @param con the {@link AlphabetContainer} of the PMM
	 * @param length the length of the PMM, i.e., the length of the sequences that can be modeled
	 * @param ess the equivalent sample size (ess) of the PMM, if 0 (zero) the model can be trained using maximum likelihood principle, otherwise it can be trained using the maximum a posteriori principle using the BDeu prior
	 * @param order the order of the PMM, i.e., the number of random variables (=positions) that might have an influence on the probability of outcome of a random variable (=position) 
	 * 
	 * @return the PMM
	 * 
	 * @throws Exception if the model can not be created correctly
	 */
	public static BayesianNetworkModel createPermutedMarkovModel( AlphabetContainer con, int length, double ess, byte order ) throws Exception {
		return createBN( con, length, ess, ModelType.PMM, order );
	}
	
	/**
	 * This method returns a Bayesian network model (BN) with user-specified order. Bayesian network determine a directed acyclic graph of the random variables (positions) and than learn the parameters of the distribution based on this network.
	 * 
	 * @param con the {@link AlphabetContainer} of the BN
	 * @param length the length of the BN, i.e., the length of the sequences that can be modeled
	 * @param ess the equivalent sample size (ess) of the BN, if 0 (zero) the model can be trained using maximum likelihood principle, otherwise it can be trained using the maximum a posteriori principle using the BDeu prior
	 * @param order the order of the BN, i.e., the number of random variables (=positions) that might have an influence on the probability of outcome of a random variable (=position) 
	 * 
	 * @return the BN
	 * 
	 * @throws Exception if the model can not be created correctly
	 */
	public static BayesianNetworkModel createBayesianNetworkModel( AlphabetContainer con, int length, double ess, byte order ) throws Exception {
		return createBN( con, length, ess, ModelType.BN, order );
	}
	
	/**
	 * This method returns a homogeneous Markov model with user-specified order.
	 * 
	 * @param con the {@link AlphabetContainer} of the model
	 * @param ess the equivalent sample size (ess) of the model, if 0 (zero) the model can be trained using maximum likelihood principle, otherwise it can be trained using the maximum a posteriori principle using the BDeu prior
	 * @param order the order of the model, i.e., the number of directly preceding random variables (=positions) that might have an influence on the probability of outcome of a random variable (=position) 
	 * 
	 * @return the homogeneous Markov model
	 * 
	 * @throws Exception if the model can not be created correctly
	 */
	public static HomogeneousMM createHomogeneousMarkovModel( AlphabetContainer con, double ess, byte order ) throws Exception {
		HomMMParameterSet ps = new HomMMParameterSet( con, ess, null, order );
		return (HomogeneousMM) getModelFromParameterSet( ps );
	}
	
	/**
	 * This method allows to create a {@link StrandModel} that allows to score binding sites on both strand of DNA.
	 * 
	 * @param model the internally used model
	 * 
	 * @return the {@link StrandModel}
	 * 
	 * @throws Exception if the model can not be created correctly
	 */
	public static StrandModel createStrandModel( Model model ) throws Exception {
		return new StrandModel( model, 10, 0.5, 1, new SmallDifferenceOfFunctionEvaluationsCondition(1E-6), Parameterization.LAMBDA );
	}
	
	/**
	 * This method allows to create a {@link MixtureModel} that allows to model a {@link de.jstacs.data.Sample} as a mixture of individual components.
	 * 
	 * @param hyper the hyper parameters for the components (should be identical to the ESS of the components)
	 * @param model the internally used model
	 * 
	 * @return the {@link MixtureModel}
	 * 
	 * @throws Exception if the model can not be created correctly
	 */
	public static MixtureModel createMixtureModel( double[] hyper, Model[] model ) throws Exception {
		//in most cases length can be determined by the components
		int i = 0;
		while( model[i].getLength() == 0 ) {
			i++;
		}
		return new MixtureModel( i==model.length ? 0 : model[i].getLength(), model, 10, hyper, 1, new SmallDifferenceOfFunctionEvaluationsCondition(1E-6), Parameterization.LAMBDA );
	}
	
	/**
	 * This method allows to create a &quot;zero or one occurrence per sequence&quot; (ZOOPS) model that allows to discovers binding sites in a {@link de.jstacs.data.Sample}.
	 * 
	 * @param motif the internally used model for the binding sites
	 * @param bg the internally used model for the flanking sequence
	 * @param hyper the hyper parameters for the components (should be identical to the ESS of the components)
	 * @param trainOnlyMotifModel a switch allowing to train either the motif model or both (motif and bg) models 
	 * 
	 * @return the ZOOPS model
	 * 
	 * @throws Exception if the model can not be created correctly
	 */
	public static SingleHiddenMotifMixture createZOOPS( Model motif, Model bg, double[] hyper, boolean trainOnlyMotifModel ) throws Exception {
		return new SingleHiddenMotifMixture( motif, bg, trainOnlyMotifModel, 10, hyper, null, 1, new SmallDifferenceOfFunctionEvaluationsCondition(1E-6), Parameterization.LAMBDA );
	}
}
