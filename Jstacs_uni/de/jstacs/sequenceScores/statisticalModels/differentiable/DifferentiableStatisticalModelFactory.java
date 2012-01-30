package de.jstacs.sequenceScores.statisticalModels.differentiable;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.BayesianNetworkDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.InhomogeneousMarkov;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.MixtureDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.StrandDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.StrandDiffSM.InitMethod;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.motif.ExtendedZOOPSDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.motif.UniformDurationDiffSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.ConstraintManager;

/**
 * This class allows to easily create some frequently used models.
 * It offers only one way of creating each model and set some of the parameters to default values.
 * If you like to set further models please check the constructors of the individual classes. 
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class DifferentiableStatisticalModelFactory {

	/**
	 * This method returns a position weight matrix (PWM). A PWM assumes that all positions of a sequence are statistically independent.
	 * 
	 * @param con the {@link AlphabetContainer} of the PWM
	 * @param length the length of the PWM, i.e., the length of the sequences that can be modeled
	 * @param ess the equivalent sample size (ess) of the PWM for the BDeu prior on its parameters
	 * 
	 * @return the PWM
	 * 
	 * @throws Exception if the model can not be created correctly
	 */
	public static BayesianNetworkDiffSM createPWM( AlphabetContainer con, int length, double ess ) throws Exception{
		return createInhomogeneousMarkovModel( con, length, ess, 0 );
	}
	
	/**
	 * This method returns a inhomogeneous Markov model (IMM) with user-specified order.
	 * 
	 * @param con the {@link AlphabetContainer} of the IMM
	 * @param length the length of the IMM, i.e., the length of the sequences that can be modeled
	 * @param ess the equivalent sample size (ess) of the IMM for the BDeu prior on its parameters
	 * @param order the order of the IMM, i.e., the number of directly preceding random variables (=positions) that might have an influence on the probability of outcome of a random variable (=position) 
	 * 
	 * @return the IMM
	 * 
	 * @throws Exception if the model can not be created correctly
	 */
	public static BayesianNetworkDiffSM createInhomogeneousMarkovModel( AlphabetContainer con, int length, double ess, int order ) throws Exception{
		return new BayesianNetworkDiffSM( con, length, ess, true, new InhomogeneousMarkov(order) );
	}
	
	/**
	 * This method returns a homogeneous Markov model with user-specified order.
	 * 
	 * @param con the {@link AlphabetContainer} of the model
	 * @param ess the equivalent sample size (ess) of the class of this model, used for the BDeu prior on its parameters in conjunction with <code>priorLength</code>
	 * @param order the order of the model, i.e., the number of directly preceding random variables (=positions) that might have an influence on the probability of outcome of a random variable (=position) 
	 * @param priorLength the a-priorily expected length of input sequences, is multiplied by <code>ess</code> before computing hyper-parameters
	 * 
	 * @return the homogeneous Markov model
	 */
	public static HomogeneousMMDiffSM createHomogeneousMarkovModel( AlphabetContainer con, double ess, int order, int priorLength ){
		return new HomogeneousMMDiffSM( con, order, ess, priorLength );
	}
	
	/**
	 * This method allows to create a {@link StrandDiffSM} that allows to score binding sites on both strand of DNA.
	 * The strand preferences is learned together with the supplied {@link DifferentiableStatisticalModel}.
	 * 
	 * @param model the internally used model
	 * 
	 * @return the {@link StrandDiffSM}
	 * 
	 * @throws CloneNotSupportedException if the supplied {@link DifferentiableStatisticalModel} could not be cloned
	 * @throws WrongAlphabetException if the {@link AlphabetContainer} of the supplied {@link DifferentiableStatisticalModel} is not {@link AlphabetContainer#isReverseComplementable()}
	 */
	public static StrandDiffSM createStrandModel(DifferentiableStatisticalModel model) throws CloneNotSupportedException, WrongAlphabetException{
		return new StrandDiffSM( model, 0.5, 10, true, InitMethod.INIT_BOTH_STRANDS );
	}
	
	/**
	 * This method allows to create a {@link MixtureDiffSM} that models a mixture of individual component {@link DifferentiableStatisticalModel}s.
	 * 
	 * @param models the internally used models
	 * 
	 * @return the {@link MixtureDiffSM}
	 * 
	 * @throws CloneNotSupportedException if the supplied {@link DifferentiableStatisticalModel}s could not be cloned
	 */
	public static MixtureDiffSM createMixtureModel(DifferentiableStatisticalModel[] models) throws CloneNotSupportedException{
		return new MixtureDiffSM( 10, true, models );
	}
	
	/**
	 * This method allows to create a &quot;zero or one occurrence per sequence&quot; (ZOOPS) model that allows to discover binding sites in a {@link de.jstacs.data.DataSet}.
	 * 
	 * @param length the length of the input {@link Sequence}s (only fixed length allowed)
	 * @param motif the internally used model for the binding sites
	 * @param bg the internally used model for the flanking sequence
	 * 
	 * @return the ZOOPS model
	 * 
	 * @throws Exception if the model can not be created correctly
	 */
	public static ExtendedZOOPSDiffSM createZOOPS(int length, DifferentiableStatisticalModel motif, HomogeneousDiffSM bg) throws Exception{
		return new ExtendedZOOPSDiffSM( ExtendedZOOPSDiffSM.CONTAINS_SOMETIMES_A_MOTIF, length, 10, true, bg, motif, new UniformDurationDiffSM( 0, length-motif.getLength()+1 ), true );
	}
	
	/**
	 * This method allows to create a {@link MarkovRandomFieldDiffSM} of the specified length and with the given constraint type.
	 * 
	 * @param con the {@link AlphabetContainer} of the {@link MarkovRandomFieldDiffSM}
	 * @param length the length of the {@link Sequence}s the {@link MarkovRandomFieldDiffSM} can handle
	 * @param constraintType the constraint type, see {@link ConstraintManager#extract(int, String)}
	 *
	 * @return the new {@link MarkovRandomFieldDiffSM}
	 * @see ConstraintManager#extract(int, String)
	 */
	public static MarkovRandomFieldDiffSM createMarkovRandomField(AlphabetContainer con, int length, String constraintType){
		return new MarkovRandomFieldDiffSM( con, length, constraintType );
	}
	
}
