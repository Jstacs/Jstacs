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

package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import de.jstacs.NonParsableException;
import de.jstacs.data.DataSet;
import de.jstacs.io.FileManager;
import de.jstacs.io.XMLParser;
import de.jstacs.sampling.GibbsSamplingModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.FSDAGTrainSMForGibbsSamplingParameterSet;

/**
 * This is the class for a fixed structure directed acyclic graphical model (see
 * {@link FSDAGTrainSM}) that can be used in a Gibbs sampling.
 * 
 * @author Berit Haldemann, Jens Keilwagen
 * 
 * @see de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM
 */

public class FSDAGModelForGibbsSampling extends FSDAGTrainSM implements GibbsSamplingModel {

	/**
	 * The files for saving the parameters during the sampling.
	 */
	protected File[] paramsFile;

	/**
	 * The counter for the sampling steps of each sampling.
	 */
	protected int[] counter;

	/**
	 * The index of the current sampling.
	 */
	protected int samplingIndex;

	/**
	 * The writer for the <code>paramsFile</code> in a sampling.
	 */
	protected BufferedWriter writer;

	/**
	 * The reader for the <code>paramsFile</code> after a sampling.
	 */
	protected BufferedReader reader;

	/**
	 * The default constructor.
	 * 
	 * @param params
	 *            the parameter set
	 * 
	 * @throws CloneNotSupportedException
	 *             if the parameter set could not be cloned
	 * @throws IllegalArgumentException
	 *             if the parameter set is not instantiated
	 * @throws NonParsableException
	 *             if the parameter set is not parsable
	 */
	public FSDAGModelForGibbsSampling( FSDAGTrainSMForGibbsSamplingParameterSet params ) throws CloneNotSupportedException,
																						IllegalArgumentException, NonParsableException {
		super( params );
		paramsFile = null;
		counter = null;
		reader = null;
		writer = null;
	}

	/**
	 * This is the constructor for the {@link de.jstacs.Storable} interface.
	 * Creates a new {@link FSDAGModelForGibbsSampling} out of its XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} could not be parsed.
	 */
	public FSDAGModelForGibbsSampling( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/**
	 * In this method the <code>reader</code> is set to <code>null</code> and
	 * the <code>paramsFile</code> is cloned.
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.DAGTrainSM#clone()
	 * @see java.lang.Object#clone()
	 */
	@Override
	public FSDAGModelForGibbsSampling clone() throws CloneNotSupportedException {
		FSDAGModelForGibbsSampling clone = (FSDAGModelForGibbsSampling)super.clone();
		if( writer != null ) {
			throw new CloneNotSupportedException( "sampling was not stopped before" );
		} else {
			clone.writer = null;
		}
		clone.reader = null;

		if( paramsFile != null ) {
			try {
				clone.paramsFile = new File[paramsFile.length];
				clone.counter = new int[paramsFile.length];
				for( int i = 0; i < paramsFile.length; i++ ) {
					if( paramsFile[i] != null ) {
						clone.paramsFile[i] = File.createTempFile( "fsdag-", ".dat", null );
						FileManager.copy( paramsFile[i].getAbsolutePath(), clone.paramsFile[i].getAbsolutePath() );
						clone.counter[i] = counter[i];
					}
				}
			} catch ( IOException e ) {
				CloneNotSupportedException c = new CloneNotSupportedException( e.getMessage() );
				c.setStackTrace( e.getStackTrace() );
				throw c;
			}

		}
		return clone;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sampling.SamplingComponent#parseNextParameterSet()
	 */
	public boolean parseNextParameterSet() {
		if( writer != null ) {
			return false;
		}
		String str = null;
		try {
			str = reader.readLine();
		} catch ( IOException e ) {} finally {
			if( str == null ) {
				return false;
			}
		}

		parse( str );
		return true;
	}

	private void parse( String str ) {
		String[] strArray = str.split( "\t" );
		for( int k = 1, l = 0; l < length; l++ ) {
			constraints[l].setFreqs( strArray, k );
			k += constraints[l].getNumberOfSpecificConstraints();
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sampling.SamplingComponent#parseParameterSet(int, int)
	 */
	public boolean parseParameterSet( int sampling, int n ) throws NumberFormatException, IOException {
		String str;
		if( reader != null ) {
			reader.close();
		}
		reader = new BufferedReader( new FileReader( paramsFile[sampling] ) );
		while( ( str = ( reader.readLine() ) ) != null ) {
			if( Integer.parseInt( str.substring( 0, str.indexOf( "\t" ) ) ) == n ) {
				parse( str );
				return true;
			}
		}
		return false;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sampling.SamplingComponent#initModelForSampling(int)
	 */
	public void initForSampling( int starts ) throws IOException {
		if( paramsFile != null && paramsFile.length == starts ) {
			FileOutputStream o;
			for( int i = 0; i < starts; i++ ) {
				if( paramsFile[i] != null ) {
					o = new FileOutputStream( paramsFile[i] );
					o.close();
				}
				counter[i] = 0;
			}
		} else {
			deleteParameterFiles();
			paramsFile = new File[starts];
			counter = new int[starts];
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sampling.SamplingComponent#extendSampling(int, boolean)
	 */
	public void extendSampling( int sampling, boolean extend ) throws IOException {
		if( paramsFile[sampling] == null ) {
			paramsFile[sampling] = File.createTempFile( "fsdag-", ".dat", null );
		} else {
			if( extend ) {
				parseParameterSet( sampling, counter[sampling] - 1 );
				reader.close();
				reader = null;
			} else {
				counter[sampling] = 0;
			}
		}
		writer = new BufferedWriter( new FileWriter( paramsFile[sampling], extend ) );
		samplingIndex = sampling;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.DAGTrainSM#drawParameters(de.jstacs.data.DataSet, double[])
	 */
	@Override
	public void drawParameters( DataSet data, double[] weights ) throws Exception {
		super.drawParameters( data, weights );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sampling.SamplingComponent#samplingStopped()
	 */
	public void samplingStopped() throws IOException {
		if( writer != null ) {
			writer.close();
			writer = null;
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.DAGTrainSM#getFurtherModelInfos()
	 */
	@Override
	protected StringBuffer getFurtherModelInfos() {
		if( writer != null ) {
			throw new RuntimeException( "could not parse the model to XML while sampling" );
		}
		StringBuffer xml = super.getFurtherModelInfos();
		if( xml == null ) {
			xml = new StringBuffer( 1000 );
		}
		XMLParser.appendObjectWithTags( xml, paramsFile != null, "hasParameters" );
		if( paramsFile != null ) {
			try {
				XMLParser.appendObjectWithTags( xml, counter, "counter" );
				String content;
				for( int i = 0; i < paramsFile.length; i++ ) {
					if( paramsFile[i] != null ) {
						content = FileManager.readFile( paramsFile[i] ).toString();
					} else {
						content = "";
					}
					XMLParser.appendObjectWithTagsAndAttributes( xml, content, "fileContent", "pos=\"" + i + "\"" );
				}
			} catch ( IOException e ) {
				RuntimeException r = new RuntimeException( e.getMessage() );
				r.setStackTrace( e.getStackTrace() );
				throw r;
			}
		}
		return xml;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.DAGTrainSM#setFurtherModelInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void setFurtherModelInfos( StringBuffer xml ) throws NonParsableException {
		super.setFurtherModelInfos( xml );
		if( XMLParser.extractObjectForTags( xml, "hasParameters", boolean.class ) ) {
			counter = XMLParser.extractObjectForTags( xml, "counter", int[].class );
			paramsFile = new File[counter.length];
			try {
				String content;
				Map<String,String> filter = new TreeMap<String, String>();
				for( int i = 0; i < paramsFile.length; i++ ) {
					filter.clear();
					filter.put( "pos", ""+i );
					content = XMLParser.extractObjectAndAttributesForTags( xml, "fileContent", null, filter, String.class );
					if( !content.equalsIgnoreCase( "" ) ) {
						paramsFile[i] = File.createTempFile( "pi-", ".dat", null );
						FileManager.writeFile( paramsFile[i], new StringBuffer( content ) );
					}
				}
			} catch ( IOException e ) {
				NonParsableException n = new NonParsableException( e.getMessage() );
				n.setStackTrace( e.getStackTrace() );
				throw n;
			}
		} else {
			counter = null;
			paramsFile = null;
		}
		writer = null;
		reader = null;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sampling.SamplingComponent#isInSamplingMode()
	 */
	public boolean isInSamplingMode() {
		return writer != null;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.FSDAGTrainSM#drawParameters(de.jstacs.data.DataSet, double[], int[][])
	 */
	@Override
	public void drawParameters( DataSet data, double[] weights, int[][] graph ) throws Exception {
		if( isInSamplingMode() ) {
			throw new RuntimeException( "could not change the structure while sampling" );
		}
		super.drawParameters( data, weights, graph );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.FSDAGTrainSM#train(de.jstacs.data.DataSet, double[])
	 */
	@Override
	public void train( DataSet data, double[] weights ) throws Exception {
		if( isInSamplingMode() ) {
			throw new RuntimeException( "could not train the model while sampling" );
		}
		super.train( data, weights );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.FSDAGTrainSM#train(de.jstacs.data.DataSet, double[], int[][])
	 */
	@Override
	public void train( DataSet data, double[] weights, int[][] graph ) throws Exception {
		if( isInSamplingMode() ) {
			throw new RuntimeException( "could not train the model while sampling" );
		}
		super.train( data, weights, graph );
	}

	/* 
	 * @see java.lang.Object#finalize()
	 */
	@Override
	protected void finalize() throws Throwable {
		if( writer != null ) {
			writer.close();
		}
		if( reader != null ) {
			reader.close();
		}
		deleteParameterFiles();
		super.finalize();
	}

	private void deleteParameterFiles() {
		if( paramsFile != null ) {
			for( int i = 0; i < paramsFile.length; i++ ) {
				if( paramsFile[i] != null ) {
					paramsFile[i].delete();
				}
			}
		}
	}

    public void acceptParameters() throws IOException {
		writer.write("" + (counter[samplingIndex]++));
		for (int n, i, l = 0; l < length; l++) {
			n = constraints[l].getNumberOfSpecificConstraints();
			for (i = 0; i < n; i++) {
				writer.write("\t" + constraints[l].getFreq(i));
			}
		}
		writer.write("\n");
		writer.flush();
	}
}
