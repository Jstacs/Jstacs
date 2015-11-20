package de.jstacs.sequenceScores.statisticalModels.trainable;

import java.text.NumberFormat;

import de.jstacs.clustering.hierachical.PWMSupplier;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.utils.Normalisation;


public class PFMWrapperTrainSM extends AbstractTrainableStatisticalModel implements PWMSupplier {

	private double[][] logPWM;
	private double[][] pfm;
	private String name;
	
	public PFMWrapperTrainSM( AlphabetContainer alphabets, String name, double[][] pfm, double ess ) throws CloneNotSupportedException {
		super( alphabets, pfm.length );
		//TODO check
		this.pfm = (double[][])ArrayHandler.clone(pfm);
		logPWM = new double[pfm.length][];
		for(int i=0;i<pfm.length;i++){
			logPWM[i] = new double[pfm[i].length];
			for(int j=0;j<logPWM[i].length;j++){
				logPWM[i][j] = Math.log( pfm[i][j] + ess / pfm[i].length );
			}
			double norm = Normalisation.getLogSum( logPWM[i] );
			for(int j=0;j<logPWM[i].length;j++){
				logPWM[i][j] -= norm;
			}
		}
		this.name = name;
	}

	public PFMWrapperTrainSM( StringBuffer stringBuff ) throws NonParsableException {
		super( stringBuff );
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer sb = new StringBuffer();
		XMLParser.appendObjectWithTags( sb, alphabets, "alphabet" );
		XMLParser.appendObjectWithTags( sb, logPWM, "logPWM" );
		XMLParser.appendObjectWithTags(sb, pfm, "pfm");
		XMLParser.appendObjectWithTags( sb, name, "name" );
		XMLParser.addTags( sb, getClass().getSimpleName() );
		return sb;
	}

	@Override
	public void train( DataSet data, double[] weights ) throws Exception {
		
	}

	@Override
	public double getLogProbFor( Sequence sequence, int startpos, int endpos ) throws Exception {
		double prob = 0;
		for(int i=startpos;i<=endpos;i++){
			prob += logPWM[i-startpos][sequence.discreteVal( i )];
		}
		return prob;
	}

	@Override
	public double getLogPriorTerm() throws Exception {
		return 0;
	}

	@Override
	public String getInstanceName() {
		return name == null ? "PFM Wrapper" : name;
	}
	
	@Override
	public String getName() {
		return name;
	}

	@Override
	public NumericalResultSet getNumericalCharacteristics() throws Exception {
		return null;
	}

	@Override
	public boolean isInitialized() {
		return true;
	}

	
	@Override
	public String toString( NumberFormat nf ) {
		if(name != null){
			return name;
		}else{
			StringBuffer sb = new StringBuffer();
			for(int i=0;i<logPWM.length;i++){
				sb.append( nf.format( Math.exp( logPWM[i][0] ) ) );
				for(int j=1;j<logPWM[i].length;j++){
					sb.append( "\t" + nf.format( Math.exp( logPWM[i][0] ) ) );
				}
				sb.append( "\n" );
			}
			return sb.toString();
		}
	}

	@Override
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, getClass().getSimpleName() );
		alphabets = (AlphabetContainer)XMLParser.extractObjectForTags( xml, "alphabet" );
		logPWM = (double[][])XMLParser.extractObjectForTags( xml, "logPWM" );
		pfm = (double[][]) XMLParser.extractObjectForTags(xml, "pfm");
		name = (String)XMLParser.extractObjectForTags( xml, "name" );
		this.length = logPWM.length;
	}
	
	public double[][] getPFM() throws CloneNotSupportedException{
		return ArrayHandler.clone(pfm);
	}

	public double[][] getPWM() {
		double[][] pwm = new double[logPWM.length][];
		for(int i=0;i<pwm.length;i++){
			pwm[i] = new double[logPWM[i].length];
			for(int j=0;j<pwm[i].length;j++){
				pwm[i][j] = Math.exp( logPWM[i][j] );
			}
		}
		return pwm;
	}

}
