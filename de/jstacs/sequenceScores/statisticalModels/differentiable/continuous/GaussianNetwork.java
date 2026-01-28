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
package de.jstacs.sequenceScores.statisticalModels.differentiable.continuous;

import java.text.NumberFormat;
import java.util.Arrays;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.random.RandomNumberGenerator;

/**
 * Class for a fixed-structure Gaussian network, i.e., the structure of the Gaussian network needs to be provided at construction
 * and remains fixed. Conditional Gaussians at each node of the network are parameterized in terms of log-precision, mean, and dependency parameters
 * b.
 * The local likelihood at node \( i \) with parents \( pa(i) \) is defined as
 * $$ \begin{equation} P(x_i | pa(i)) = \sqrt{ \frac{\exp(\lambda)}{2 \pi} } \cdot \exp\left( -\frac{ \exp(\lambda) }{2} \left( x_i - \mu_i - \sum_{j \in pa(i)} b_{ij}\left( x_j - \mu_j \right)  \right)^2 \right) \end{equation}.&#92;$$
 * Currently, this {@link de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel} is defined only for non-Bayesian learning, since {@link #addGradientOfLogPriorTerm(double[], int)}
 * and {@link #getESS()} are not implemented reasonably.
 *   
 * @author Jan Grau
 *
 */
public class GaussianNetwork extends AbstractDifferentiableStatisticalModel {
	
	private double[] mu;
	private double[] lambda;
	private double[][] bij;
	private int[][] structure;
	private int[] boff;
	private double ess;
	
	/**
	 * Create a new Gaussian network with the given <code>structure</code>.
	 * @param structure the parents per position, provided as a two-dimensional array, 
	 * where the first dimension indexes the nodes while the second dimension holds the indexes of the parent variables, e.g.,
	 * <code>structure[1]</code> holds the parents of node 1, and <code>structure[1] = {2,3}</code> denotes that node has has parents 2 and 3.
	 * May also be empty, i.e., <code>structure[1]</code>, which results in only unconditional probabilities at this position.
	 * @throws CloneNotSupportedException if the structure array could not be cloned
	 */
	public GaussianNetwork(int[][] structure) throws CloneNotSupportedException {
		this(new AlphabetContainer(new ContinuousAlphabet()), structure );
	}
	
	/**
	 * Create a new Gaussian network with the given <code>structure</code> and an explicitly defined {@link AlphabetContainer}.
	 * @param con the alphabet
	 * @param structure the parents per position, provided as a two-dimensional array, 
	 * where the first dimension indexes the nodes while the second dimension holds the indexes of the parent variables, e.g.,
	 * <code>structure[1]</code> holds the parents of node 1, and <code>structure[1] = {2,3}</code> denotes that node has has parents 2 and 3.
	 * May also be empty, i.e., <code>structure[1]</code>, which results in only unconditional probabilities at this position.
	 * @throws CloneNotSupportedException if the structure array could not be cloned
	 */
	public GaussianNetwork(AlphabetContainer con, int[][] structure) throws CloneNotSupportedException {
		super(con, structure.length);
		this.structure = ArrayHandler.clone(structure);
		this.mu = new double[structure.length];
		this.lambda = new double[structure.length];
		bij = new double[structure.length][];
		int off = mu.length+lambda.length;
		boff = new int[bij.length];
		for(int i=0;i<structure.length;i++){
			boff[i] = off;
			bij[i] = new double[structure[i].length];
			off += bij[i].length;
		}
	}
	
	/**
	 * Create a new Gaussian network with the given <code>structure</code> and an explicitly defined {@link AlphabetContainer}.
	 * @param con the alphabet
	 * @param structure the parents per position, provided as a two-dimensional array, 
	 * where the first dimension indexes the nodes while the second dimension holds the indexes of the parent variables, e.g.,
	 * <code>structure[1]</code> holds the parents of node 1, and <code>structure[1] = {2,3}</code> denotes that node has has parents 2 and 3.
	 * May also be empty, i.e., <code>structure[1]</code>, which results in only unconditional probabilities at this position.
	 * @param ess the equivalent sample size, currently not used in a parameter prior
	 * @throws CloneNotSupportedException if the structure array could not be cloned
	 */
	public GaussianNetwork(AlphabetContainer con, int[][] structure, double ess) throws CloneNotSupportedException {
		this(con,structure);
		this.ess=ess;
	}
	
	
	public GaussianNetwork clone() throws CloneNotSupportedException{
		GaussianNetwork clone = (GaussianNetwork) super.clone();
		clone.structure = structure.clone();
		for(int i=0;i<clone.structure.length;i++){
			clone.structure[i] = structure[i].clone();
		}
		clone.mu = mu.clone();
		clone.lambda = lambda.clone();
		clone.bij = bij.clone();
		for(int i=0;i<clone.bij.length;i++){
			clone.bij[i] = bij[i].clone();
		}
		clone.boff = boff.clone();
		return clone;
	}
	
	/**
	 * Creates a Gaussian network from its XML representation.
	 * @param xml the XML representation
	 * @throws NonParsableException if <code>xml</code> could not be parsed
	 */
	public GaussianNetwork(StringBuffer xml) throws NonParsableException {
		super(xml);
	}

	@Override
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		return 0;
	}

	@Override
	public double getLogNormalizationConstant() {
		return 0;
	}

	@Override
	public double getLogPartialNormalizationConstant(int parameterIndex) throws Exception {
		return Double.NEGATIVE_INFINITY;
	}

	@Override
	public double getLogPriorTerm() {
		return 0;
	}

	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int start) throws Exception {

	}

	@Override
	public double getESS() {
		return ess;
	}

	@Override
	public void initializeFunction(int index, boolean freeParams, DataSet[] data, double[][] weights) throws Exception {
		
		initializeFunctionRandomly(freeParams);
		
		DataSet dat = data[index];
		
		Arrays.fill(mu, 0);
		Arrays.fill(lambda, 0);
		
		double n = 0;
		
		for(int i=0;i<dat.getNumberOfElements();i++){
			double w = weights == null || weights[index] == null ? 1 : weights[index][i];
			
			Sequence seq = dat.getElementAt(i);
			for(int j=0;j<seq.getLength();j++){
				double val = seq.continuousVal(j);
				mu[j] += val*w;
				lambda[j] += val*val*w;
			}
			n += w;
		}
		
		for(int i=0;i<mu.length;i++){
			mu[i] /= n;
			lambda[i] /= n;
			//System.out.println(i+" "+mu[i]+" "+lambda[i]+" "+n);
			
			lambda[i] = lambda[i] - mu[i]*mu[i];
			lambda[i] = -Math.log(lambda[i]);
		}
		
	}

	@Override
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
		RandomNumberGenerator rng = new RandomNumberGenerator();
		for(int i=0;i<mu.length;i++){
			mu[i] = rng.nextGaussian();
			lambda[i] = Math.log(0.1);//rng.nextGammaLog(1, 1);
			for(int j=0;j<bij[i].length;j++){
				bij[i][j] = rng.nextGaussian();
			}
		}
	}

	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList partialDer) {
		
		double val = 0.0;
				
		for(int i=0;i<structure.length;i++){
			
			double mymu = mu[i];
			for(int j=0;j<structure[i].length;j++){
				mymu += bij[i][j]*( seq.continuousVal(start+structure[i][j]) - mu[structure[i][j]] );
			}
			double temp = (seq.continuousVal(start+i) - mymu);
			double expl = Math.exp(lambda[i]);
			val += 0.5*lambda[i] - 0.5*Math.log(2.0*Math.PI) - expl/2.0*temp*temp;
			
			//lamnda_i
			indices.add(mu.length+i);
			partialDer.add( 0.5 - expl/2.0*temp*temp );
			
			//mu_i
			indices.add(i);
			partialDer.add( expl*temp );
			
			//mu_k
			for(int j=0;j<structure[i].length;j++){
				indices.add(structure[i][j]);
				partialDer.add( -expl*temp*bij[i][j] );
			}
			
			//bij
			for(int j=0;j<structure[i].length;j++){
				indices.add(boff[i]+j);
				partialDer.add( expl*temp*( seq.continuousVal(start+structure[i][j]) - mu[structure[i][j]] ) );
			}
			
		}
		
		int derEnd = partialDer.length();
		
		return val;
		
	}

	@Override
	public int getNumberOfParameters() {
		return boff[boff.length-1] + bij[bij.length-1].length;
	}

	@Override
	public double[] getCurrentParameterValues() throws Exception {
		double[] res = new double[getNumberOfParameters()];
		int off = 0;
		System.arraycopy(mu, 0, res, 0, mu.length);
		off += mu.length;
		System.arraycopy(lambda, 0, res, off, lambda.length);
		off += lambda.length;
		for(int i=0;i<bij.length;i++){
			System.arraycopy(bij[i], 0, res, off, bij[i].length);
			off += bij[i].length;
		}
		return res;
	}

	@Override
	public void setParameters(double[] params, int start) {
		System.arraycopy(params, start, mu, 0, mu.length);
		start += mu.length;
		System.arraycopy(params, start, lambda, 0, lambda.length);
		start += lambda.length;
		for(int i=0;i<bij.length;i++){
			System.arraycopy(params, start, bij[i], 0, bij[i].length);
			start += bij[i].length;
		}
	}

	@Override
	public String getInstanceName() {
		return "GaussianNetwork";
	}

	@Override
	public double getLogScoreFor(Sequence seq, int start) {
		
		double val = 0.0;
		
		for(int i=0;i<structure.length;i++){
			
			double mymu = mu[i];
			for(int j=0;j<structure[i].length;j++){
				mymu += bij[i][j]*( seq.continuousVal(start+structure[i][j]) - mu[structure[i][j]] );
			}
			double temp = (seq.continuousVal(start+i) - mymu);
			val += 0.5*lambda[i] - 0.5*Math.log(2.0*Math.PI) - Math.exp(lambda[i])/2.0*temp*temp;
		}
		//System.out.println("g\t"+val);
		return val;
	}

	@Override
	public boolean isInitialized() {
		return true;
	}

	@Override
	public String toString(NumberFormat nf) {
		StringBuffer sb = new StringBuffer();
		//sb.append("mus: "+Arrays.toString(mu)+"\n");
		//sb.append("lambdas: "+Arrays.toString(lambda)+"\n");
		for(int i=0;i<mu.length;i++){
			sb.append(i+": dnorm("+mu[i]+","+(1.0/Math.exp(lambda[i]))+")\n");
		}
		for(int i=0;i<bij.length;i++){
			sb.append("b_"+i+": "+Arrays.toString(bij[i])+"\n");
		}
		return sb.toString();
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, bij, "bij");
		XMLParser.appendObjectWithTags(xml, boff, "boff");
		XMLParser.appendObjectWithTags(xml, lambda, "lambda");
		XMLParser.appendObjectWithTags(xml, mu, "mu");
		XMLParser.appendObjectWithTags(xml, structure, "structure");
		XMLParser.appendObjectWithTags(xml, ess, "ess");
		XMLParser.addTags(xml, "GaussNet");
		return xml;
	}

	@Override
	protected void fromXML(StringBuffer xml) throws NonParsableException {
		xml = XMLParser.extractForTag(xml, "GaussNet");
		bij = (double[][]) XMLParser.extractObjectForTags(xml, "bij");
		boff = (int[]) XMLParser.extractObjectForTags(xml, "boff");
		lambda = (double[]) XMLParser.extractObjectForTags(xml, "lambda");
		mu = (double[]) XMLParser.extractObjectForTags(xml, "mu");
		structure = (int[][]) XMLParser.extractObjectForTags(xml, "structure");
		try {
			ess=(Double) XMLParser.extractObjectForTags(xml, "ess") ;
		} catch (NonParsableException ex) {
			ess = 0.0;
		}
		alphabets = new AlphabetContainer(new ContinuousAlphabet());
		length = mu.length;
	}

	
}
