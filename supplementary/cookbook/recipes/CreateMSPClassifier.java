package supplementary.cookbook.recipes;

import de.jstacs.classifiers.AbstractClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.msp.MSPClassifier;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModelFactory;


public class CreateMSPClassifier {

	/**
	 * @param args
	 */
	public static void main( String[] args ) throws Exception {
		AlphabetContainer con = DNAAlphabetContainer.SINGLETON;
		DifferentiableStatisticalModel pwm = DifferentiableStatisticalModelFactory.createPWM(con, 10, 4);
		GenDisMixClassifierParameterSet pars = new GenDisMixClassifierParameterSet(con,10,(byte)10,1E-9,1E-10,1,false,KindOfParameter.PLUGIN,true,1);
		AbstractClassifier cl = new MSPClassifier( pars, pwm, pwm );
	}
}