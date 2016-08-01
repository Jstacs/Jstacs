package de.jstacs.sequenceScores.statisticalModels.differentiable.continuous;

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Random;
import umontreal.iro.lecuyer.util.Num;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.ToolBox;


public class DirichletDiffSM extends AbstractDifferentiableStatisticalModel {

	private static Random r = new Random();
	
	private double logNorm;
	private double[] partNorm;
	private double ess;
	private double[] pars;
	private double[] expPars;
	private double[] chi;

	private boolean isInitialized;
	private int numStarts;
	
	
	/**
	 * @param alphabets
	 * @param length
	 * @throws IllegalArgumentException
	 */
	public DirichletDiffSM( AlphabetContainer alphabets, int length, double[] chi, double ess, int numStarts ) throws IllegalArgumentException {
		super( alphabets, length );
		if(chi.length != length+1){
			throw new IllegalArgumentException();
		}
		this.pars = new double[length+1];
		this.ess = ess;
		this.expPars = new double[length+1];
		Arrays.fill( expPars, 1 );
		this.chi = chi.clone();
		this.numStarts = numStarts;
	}
	
	public DirichletDiffSM clone() throws CloneNotSupportedException{
		DirichletDiffSM clone = (DirichletDiffSM)super.clone();
		clone.chi = chi.clone();
		clone.expPars = expPars.clone();
		clone.pars = pars.clone();
		if(partNorm != null){
			clone.partNorm = partNorm.clone();
		}
		return clone;
	}

	/**
	 * @param xml
	 * @throws NonParsableException
	 */
	public DirichletDiffSM( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public int getSizeOfEventSpaceForRandomVariablesOfParameter( int index ) {
		return 1;
	}
	

	@Override
	public int getNumberOfRecommendedStarts() {
		return numStarts;
	}

	@Override
	public double getLogNormalizationConstant() {
		return 0;
	}

	@Override
	public double getLogPartialNormalizationConstant( int parameterIndex ) throws Exception {
		return Double.NEGATIVE_INFINITY;
	}

	@Override
	public double getESS() {
		return ess;
	}

	@Override
	public double getLogPriorTerm() {
		double prior = logNorm*ess;
		for(int i=0;i<pars.length;i++){
			prior += ess*chi[i]*expPars[i] + pars[i];
		}
		return prior;
	}

	@Override
	public void addGradientOfLogPriorTerm( double[] grad, int start ) throws Exception {
		for(int i=0;i<pars.length;i++){
			grad[i+start] += ess*( partNorm[i] + chi[i]*expPars[i] ) + 1;
		}
	}

	@Override
	public void initializeFunction( int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception {
		
		double[] sums = new double[length+1];
		for(int i=0;i<data[index].getNumberOfElements();i++){
			double rest = 1.0;
			for(int j=0;j<sums.length;j++){
				double temp = 0;
				if(j<sums.length-1){
					temp = data[index].getElementAt( i ).continuousVal( j );
					rest -= temp;
					temp *= weights[index][i];
				}else{
					temp = rest*weights[index][i];
				}
				sums[j] += temp;
			}
			
		}
		
		Normalisation.sumNormalisation( sums );
		
		for(int i=0;i<sums.length;i++){
			sums[i] = Math.log( sums[i]*2+2 );
		}
		
		setParameters( sums, 0 );
		
		isInitialized = true;
		
	}

	@Override
	public void initializeFunctionRandomly( boolean freeParams ) throws Exception {
		for(int i=0;i<pars.length;i++){
			pars[i] = 1.0+Math.abs(r.nextGaussian());
		}
		precompute();
		isInitialized = true;
	}

	@Override
	public String getInstanceName() {
		return getClass().getSimpleName();
	}

	@Override
	public double getLogScoreFor( Sequence seq, int start ) {
		double score = logNorm;
		double rest = 1.0;
		for(int i=0;i<pars.length-1;i++){
		//	System.out.println(seq.continuousVal( start+i ));
			score += Math.log( seq.continuousVal( start+i ) )*(expPars[i]-1);
			rest -= seq.continuousVal( start+i );
		}
		score += Math.log( rest )*(expPars[expPars.length-1]-1);
		//System.out.println(rest);
		//System.out.println("######################");
		
		if(Double.isInfinite( score )){
			System.out.println(Arrays.toString( pars )+" "+Arrays.toString( expPars )+" "+seq);
		}
		return score;
	}

	@Override
	public double getLogScoreAndPartialDerivation( Sequence seq, int start, IntList indices, DoubleList partialDer ) {
		double score = logNorm;
		double rest = 1.0;
		for(int i=0;i<pars.length;i++){
			double temp = 0;
			if(i<pars.length-1){
				temp = Math.log( seq.continuousVal( start+i ) );
				rest -= seq.continuousVal( start+i );
			}else{
				temp = Math.log( rest );
			}

			score += temp * (expPars[i]-1);
			indices.add( i );
			partialDer.add( partNorm[i] + temp*expPars[i] );
		}
		return score;
	}

	@Override
	public int getNumberOfParameters() {
		return pars.length;
	}

	@Override
	public double[] getCurrentParameterValues() throws Exception {
		return pars.clone();
	}

	@Override
	public void setParameters( double[] params, int start ) {
		for(int i=0;i<pars.length;i++){
			pars[i] = params[i+start];
		}
		precompute();
	}

	private void precompute() {
		try{
		logNorm = 0.0;
		if(partNorm == null){
			partNorm = new double[pars.length];
		}
		
		for(int i=0;i<pars.length;i++){
			expPars[i] = Math.exp( pars[i] );
		}
		
		logNorm += Num.lnGamma( ToolBox.sum( expPars ) );
		
		for(int i=0;i<pars.length;i++){
			logNorm -= Num.lnGamma( expPars[i] );
			partNorm[i] = expPars[i]*Num.digamma( ToolBox.sum( expPars ) ) - expPars[i]*Num.digamma( expPars[i] );
		}
		}catch(Exception e){
			System.out.println(Arrays.toString( pars ));
			System.out.println(Arrays.toString( expPars ));
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public boolean isInitialized() {
		return isInitialized;
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, alphabets, "alphabets" );
		XMLParser.appendObjectWithTags( xml, length, "length" );
		XMLParser.appendObjectWithTags( xml, ess, "ess" );
		XMLParser.appendObjectWithTags( xml, pars, "pars" );
		XMLParser.appendObjectWithTags( xml, chi, "chi" );
		XMLParser.appendObjectWithTags( xml, isInitialized, "isInitialized" );
		XMLParser.appendObjectWithTags( xml, numStarts, "numStarts" );
		XMLParser.addTags( xml, "Dirichlet" );
		return xml;
	}

	@Override
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, "Dirichlet" );
		alphabets = (AlphabetContainer)XMLParser.extractObjectForTags( xml, "alphabets" );
		length = XMLParser.extractObjectForTags( xml, "length", int.class );
		ess = XMLParser.extractObjectForTags( xml, "ess", double.class );
		pars = (double[])XMLParser.extractObjectForTags( xml, "pars" );
		chi = (double[])XMLParser.extractObjectForTags( xml, "chi" );
		isInitialized = XMLParser.extractObjectForTags( xml, "isInitialized", boolean.class );
		numStarts = XMLParser.extractObjectForTags( xml, "numStarts", int.class );
		expPars = new double[pars.length];
		
		precompute();
	}
	
	public String toString(NumberFormat nf){
		return Arrays.toString( expPars );
	}

}
