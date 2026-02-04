package de.jstacs.sequenceScores.statisticalModels.differentiable;

import java.text.NumberFormat;
import java.util.Arrays;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;

/**
 * SImilar to a {@link CyclicMarkovModelDiffSM} but without hidden variable.
 * 
 * @author Jens Keilwagen
 */
public class PeriodicHomogeneousModel extends AbstractDifferentiableStatisticalModel {

	private AlphabetContainer con;
	private int a, order, np, startPhase;
	private double[][][] params;
	private double[][][] probs;
	private double[][] logNorm;
	private int[][][] stat;
	
	private static final String XML_TAG = "PHM";
	
	public PeriodicHomogeneousModel(AlphabetContainer con, int l, int period, int order ) {
		super( con, l );
		if( !con.isSimple() ) {
			throw new RuntimeException();
		}
		this.con = con;
		a = (int) con.getAlphabetLengthAt(0);
		this.order = order;
		params = new double[period][(int) Math.pow(a, order)][a];
		prepare();
	}
	
	public PeriodicHomogeneousModel( StringBuffer xml ) throws NonParsableException {
		super(xml);
	}
	
	@Override
	protected void fromXML(StringBuffer xml) throws NonParsableException {
		xml = XMLParser.extractForTag(xml, XML_TAG);
		con = (AlphabetContainer) XMLParser.extractObjectForTags(xml, "AlphabetContainer");
		a = con.getAlphabetIndexForPosition(0);
		params  = (double[][][])  XMLParser.extractObjectForTags(xml, "parameter");
		prepare();
	}
	
	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, con, "AlphabetContainer");
		XMLParser.appendObjectWithTags(xml, params, "parameter");
		XMLParser.addTags(xml, XML_TAG);
		return xml;
	}
	
	private void prepare() {
		np = params[0].length*a;
		probs = new double[params.length][params[0].length][a];
		logNorm = new double[params.length][params[0].length];
		stat = new int[params.length][params[0].length][a];
		compute();
	}
	
	private void compute() {
		for( int i = 0; i < params.length; i++ ) {
			for( int j = 0; j < params[i].length; j++ ) {
				logNorm[i][j] = Normalisation.logSumNormalisation(params[i][j],0,a,probs[i][j],0);
			}
		}
	}
	
	public PeriodicHomogeneousModel clone() throws CloneNotSupportedException {
		PeriodicHomogeneousModel clone = (PeriodicHomogeneousModel) super.clone();
		clone.params = new double[params.length][params[0].length][];
		for( int i = 0; i < params.length; i++ ) {
			for( int j = 0; j < params[i].length; j++ ) {
				clone.params[i][j] = params[i][j].clone();
			}
		}
		clone.prepare();
		return clone;
	}

	@Override
	public void initializeFunctionRandomly( boolean freeParams ) {
		//TODO
		double x = 0;//-Math.log(a);
		for( int i = 0; i < params.length; i++ ) {
			for( int j = 0; j < params[i].length; j++ ) {
				Arrays.fill( params[i][j], x );
			}
		}
		compute();
	}

	@Override
	public void initializeFunction(int index, boolean freeParams, DataSet[] data, double[][] weights) throws Exception {
		if( length>1 || params.length>1 ) {
			initializeFunctionRandomly(freeParams);
		} else {
			double[][][] counts = new double[params.length][params[0].length][a];
			for( int i = 0; i < counts.length; i++ ) {
				for( int j = 0; j < counts[0].length; j++ ) {
					Arrays.fill(counts[i][j], 1);
				}
			}
			for( int n = 0; n < data[index].getNumberOfElements(); n++ ) {
				Sequence seq = data[index].getElementAt(n);
				int idx = 0, p = startPhase;
				for( int i = order; i > 0; i-- ) {
					int d = order-i>=0 ? seq.discreteVal(order-i) : 0;
					idx = (idx*a+d);
				}
				for( int i = order; i < seq.getLength(); i++ ) {
					int d = seq.discreteVal(i);
					counts[p][idx][d]+=weights[index][n];
					idx = (idx*a+d) % params[0].length;
				}
			}
			for( int i = 0; i < counts.length; i++ ) {
				for( int j = 0; j < counts[0].length; j++ ) {
					double sum = 0;
					for( int k = 0; k < a; k++ ) {
						sum += counts[i][j][k];
					}
					if( sum>a ) {
						for( int k = 0; k < a; k++ ) {
							params[i][j][k] = Math.log( counts[i][j][k] / sum );
						}
					} else {
						System.arraycopy(params[i][0], 0, params[i][j],  0,  a );
					}
				}
			}
			compute();
		}
	}
	
	public boolean setStartPhase( int phase ) {
		if( params.length>1 ) {
			if( phase < 0 || phase >= params.length ) {
				throw new IllegalArgumentException(phase +" is not in [0, " + (params.length-1) +"]" ); 
			}
			startPhase = phase;
			return true;
		}
		return false;
	}
	
	@Override
	public double getLogScoreFor(Sequence seq, int start) {
		return getLogScoreAndPartialDerivation(seq, start, null, null);
	}

	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList partialDer) {
		int p = startPhase;
		
		//initial context
		int idx = 0;
		for( int i = order; i > 0; i-- ) {
			int d = start-i>=0 ? seq.discreteVal(start - i) : 0;
			idx = (idx*a+d);
		}
		//computation
		for( int i = 0; i < length; i++ ) {
			int d = seq.discreteVal(start + i);
			stat[p][idx][d]++;
		
			idx = (idx*a+d) % params[0].length;
			p++;
			if( p==params.length ) p = 0;
		}
		
		double res = 0;
		for( int i = 0; i < stat.length; i++ ) {
			for( int j = 0; j < stat[i].length; j++ ) {
				int s = 0;
				for( int d = 0; d < a; d++ ) {
					s += stat[i][j][d];
				}
				if( s>0 ) {
					for( int d = 0; d < a; d++ ) {
						res += stat[i][j][d]*params[i][j][d];
						if( partialDer!= null ) {
							indices.add( i*np + j*a + d );
							partialDer.add( stat[i][j][d] - s*probs[i][j][d] );
						}
						stat[i][j][d]=0;
					}
					res -= s*logNorm[i][j];
				}
			}
		}
		return res;
	}

	@Override
	public String toString(NumberFormat nf) {
		StringBuffer sb = new StringBuffer();
		DiscreteAlphabet abc = (DiscreteAlphabet) con.getAlphabetAt(0);
		for( int i = 0; i < params.length; i++ ) {
			for( int d = 0; d < a; d++ ) {
				sb.append( abc.getSymbolAt(d) + "\t" );
			}
			sb.append("\n");
			for( int j = 0; j < params[i].length; j++ ) {
				for( int d = 0; d < a; d++ ) {
					sb.append( nf.format( probs[i][j][d] ) + "\t" );
				}
				sb.append("\n");
			}
			sb.append("\n");
		}
		return sb.toString();
	}

	@Override
	public double[] getCurrentParameterValues() throws Exception {
		double[] p = new double[params.length*params[0].length*a];
		int offset=0;
		for( int i = 0; i < params.length; i++ ) {
			for( int j = 0; j < params[0].length; j++ ) {
				System.arraycopy( params[i][j], 0, p, offset, a );
				offset += a;
			}
		}
		return p;
	}

	@Override
	public void setParameters(double[] params, int start) {
		for( int i = 0; i < this.params.length; i++ ) {
			for( int j = 0; j < this.params[i].length; j++ ) {
				System.arraycopy( params, start, this.params[i][j], 0, a );
				start += a;
			}
		}
		compute();
	}

	@Override
	public String getInstanceName() {
		return "PHM " + params.length + " " + (Math.log(params[0].length)/Math.log(a));
	}

	@Override
	public boolean isInitialized() {
		return true;
	}

	@Override
	public int getNumberOfParameters() {
		return params.length*np;
	}
	

	@Override
	public boolean isNormalized() {
		return false;
	}
	
	//XXX not implemented
	@Override
	public double getLogPriorTerm() {
		throw new RuntimeException();
	}
		
	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int offset) {
		throw new RuntimeException();
	}

	@Override
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		throw new RuntimeException();
	}

	@Override
	public double getLogNormalizationConstant() {
		throw new RuntimeException();
	}

	@Override
	public double getLogPartialNormalizationConstant(int parameterIndex) throws Exception {
		throw new RuntimeException();
	}

	@Override
	public double getESS() {
		return 0;
	}
	
	public byte getMaximalMarkovOrder() {
		return (byte) order;
	}
}