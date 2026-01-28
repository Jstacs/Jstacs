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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package projects.dream2016.mix;

import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;


/**
 * 
 * @author Jens Keilwagen
 */
public interface OptimizableClassifier extends Cloneable {
	
	
	public OptimizableClassifier clone() throws CloneNotSupportedException;
	
	/**
	 * 
	 * @param classIndex the class index
	 * @param seq the sequence
	 * 
	 * @return \( \log P(classIndex | seq) \)
	 */
	double getLogProb( int classIndex, Sequence seq ) throws EvaluationException;
	
	/**
	 * 
	 * @param classIndex the class index
	 * @param seq the sequence
	 * @param indices the indices of the parameter of the partial derivation
	 * @param partialDer the partial derivation; \( \frac{\partial\log P(classIndex | seq)}{\partial \lambda_i} \)
	 * 
	 * @return \( \log P(classIndex | seq) \)
	 */
	double getLogProbAndPartialDerivations( int classIndex, Sequence seq, IntList indices, DoubleList partialDer );
	
	double getLogPriorTerm() throws DimensionException, EvaluationException;
	
	void addGradient( double[] grad, int start ) throws EvaluationException;
	
	void setParameters( double[] params, int start ) throws Exception;
	
	double[] getCurrentParameterValues( KindOfParameter kind ) throws Exception;
	
	void initialize( DataSet[] data, double[][] weights ) throws Exception;
	
	void initializeRandomly() throws Exception;
	
	/**
	 * Returns the number of parameters.
	 * 
	 * @return 
	 */
	int getNumberOfParameters();
	
	void reset() throws Exception;
	
	//void getLogProbs( Sequence seq, int start, double[] logProbs );
}