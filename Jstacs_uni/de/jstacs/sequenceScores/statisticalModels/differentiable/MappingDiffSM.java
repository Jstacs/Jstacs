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

package de.jstacs.sequenceScores.statisticalModels.differentiable;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DiscreteAlphabetMapping;
import de.jstacs.data.sequences.MappedDiscreteSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.SimpleDiscreteSequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.motifDiscovery.MotifDiscoverer;
import de.jstacs.motifDiscovery.Mutable;
import de.jstacs.motifDiscovery.MutableMotifDiscoverer;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This class implements a {@link DifferentiableStatisticalModel} that works on
 * mapped {@link Sequence}s. For instance this can be useful for protein
 * sequences to reduce the alphabet of size 20 to a smaller alphabet using for
 * instance some chemical properties of the amino acids.
 * 
 * <b>Be careful with references to {@link Sequence}s in the internal
 * {@link DifferentiableStatisticalModel}, since the {@link Sequence}s might be
 * unexpectedly mutable.</b>
 * 
 * @author Jens Keilwagen
 * 
 * @see MappedDiscreteSequence
 * @see de.jstacs.data.alphabets.DiscreteAlphabetMapping
 */
public final class MappingDiffSM extends AbstractDifferentiableStatisticalModel implements MutableMotifDiscoverer, Mutable {

	private DifferentiableStatisticalModel nsf;

	private MutableMappedDiscreteSequence mappedSeq;

	/**
	 * The main constructor creating a {@link MappingDiffSM}.
	 * 
	 * @param originalAlphabetContainer
	 *            the original {@link AlphabetContainer}
	 * @param nsf
	 *            the internally used {@link DifferentiableStatisticalModel}
	 * @param mapping
	 *            the {@link de.jstacs.data.alphabets.DiscreteAlphabetMapping}s defining the
	 *            transformation from the original {@link AlphabetContainer} to
	 *            the {@link AlphabetContainer} of the
	 *            {@link DifferentiableStatisticalModel} <code>nsf</code>
	 * 
	 * @throws WrongAlphabetException
	 *             if there is a problem with the mapping of the
	 *             {@link de.jstacs.data.alphabets.Alphabet}s
	 * @throws CloneNotSupportedException
	 *             if the {@link DifferentiableStatisticalModel} could not be
	 *             cloned
	 */
	public MappingDiffSM( AlphabetContainer originalAlphabetContainer, DifferentiableStatisticalModel nsf,
									DiscreteAlphabetMapping... mapping ) throws WrongAlphabetException,
																						CloneNotSupportedException {
		super( originalAlphabetContainer, nsf.getLength() );
		this.nsf = (DifferentiableStatisticalModel)nsf.clone();
		mappedSeq = new MutableMappedDiscreteSequence( originalAlphabetContainer, mapping );
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}. Creates a new
	 * {@link MappingDiffSM} out of a {@link StringBuffer}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public MappingDiffSM( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	public MappingDiffSM clone() throws CloneNotSupportedException {
		MappingDiffSM clone = (MappingDiffSM)super.clone();
		clone.nsf = (DifferentiableStatisticalModel)nsf.clone();
		try {
			clone.mappedSeq = new MutableMappedDiscreteSequence( alphabets, mappedSeq.getTransformations() );
			clone.mappedSeq.setSequence( mappedSeq.getOriginalSequence() );
		} catch ( WrongAlphabetException wae ) {
			// does not happen
			throw new CloneNotSupportedException( wae.getMessage() );
		}
		return clone;
	}

	private static final String XML_TAG = "MappingDiffSM";

	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, alphabets, "alphabets" );
		XMLParser.appendObjectWithTags( xml, nsf, "nsf" );
		XMLParser.appendObjectWithTags( xml, mappedSeq.getTransformations(), "transformations" );
		XMLParser.addTags( xml, XML_TAG );
		return xml;
	}

	@Override
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, XML_TAG );
		alphabets = (AlphabetContainer) XMLParser.extractObjectForTags( xml, "alphabets" );
		nsf = (DifferentiableStatisticalModel) XMLParser.extractObjectForTags( xml, "nsf" );
		length = nsf.getLength();
		try {
			mappedSeq = new MutableMappedDiscreteSequence( alphabets, XMLParser.extractObjectForTags( xml, "transformations", DiscreteAlphabetMapping[].class ) );
		} catch ( WrongAlphabetException e ) {

			throw new NonParsableException();
		}
	}

	public void addGradientOfLogPriorTerm( double[] grad, int start ) throws Exception {
		nsf.addGradientOfLogPriorTerm( grad, start );
	}

	public double getESS() {
		return nsf.getESS();
	}

	public double getLogNormalizationConstant() {
		return nsf.getLogNormalizationConstant();
	}

	public double getLogPartialNormalizationConstant( int parameterIndex ) throws Exception {
		return nsf.getLogPartialNormalizationConstant( parameterIndex );
	}

	public double getLogPriorTerm() {
		return nsf.getLogPriorTerm();
	}

	public int getSizeOfEventSpaceForRandomVariablesOfParameter( int index ) {
		return nsf.getSizeOfEventSpaceForRandomVariablesOfParameter( index );
	}

	public double[] getCurrentParameterValues() throws Exception {
		return nsf.getCurrentParameterValues();
	}

	public String getInstanceName() {
		return XML_TAG + " for " + nsf.getInstanceName();
	}

	public double getLogScoreFor( Sequence seq, int start ) {
		mappedSeq.setSequence( seq );
		return nsf.getLogScoreFor( mappedSeq, start ) - mappedSeq.getLogNumberOfPossibleOriginalSequences( start,
						length == 0 ? seq.getLength() : ( start + length ) );
	}

	public double getLogScoreAndPartialDerivation( Sequence seq, int start, IntList indices, DoubleList partialDer ) {
		mappedSeq.setSequence( seq );
		return nsf.getLogScoreAndPartialDerivation( mappedSeq, start, indices, partialDer ) - mappedSeq.getLogNumberOfPossibleOriginalSequences( start,
						length == 0 ? seq.getLength() : ( start + length ) );
	}

	public int getNumberOfParameters() {
		return nsf.getNumberOfParameters();
	}

	public void initializeFunction( int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception {
		DataSet[] mappedData = getMappedData( data );
		nsf.initializeFunction( index, freeParams, mappedData, weights );
		for( int i = 0; i < data.length; i++ ) {
			mappedData[i] = null;
		}
		mappedData = null;
		System.gc();
	}
	
	private DataSet[] getMappedData( DataSet... data ) throws WrongAlphabetException, EmptyDataSetException {
		DataSet[] mappedData = new DataSet[data.length];
		Sequence[] seqs;
		DiscreteAlphabetMapping[] transformation = mappedSeq.getTransformations();
		for( int i = 0; i < data.length; i++ ) {
			if( data[i] != null ) {
				seqs = new Sequence[data[i].getNumberOfElements()];
				for( int n = 0; n < seqs.length; n++ ) {
					seqs[n] = new MappedDiscreteSequence( (SimpleDiscreteSequence) data[i].getElementAt( n ), transformation );
				}
				mappedData[i] = new DataSet( "mapped: " + data[i].getAnnotation(), seqs );
			}
		}
		return mappedData;
	}

	public void initializeFunctionRandomly( boolean freeParams ) throws Exception {
		nsf.initializeFunctionRandomly( freeParams );
	}

	public boolean isInitialized() {
		return nsf.isInitialized();
	}

	public void setParameters( double[] params, int start ) {
		nsf.setParameters( params, start );
	}

	public String toString() {
		return getInstanceName() + "\n" + nsf.toString();
	}

	/**
	 * This method return the internal function.
	 * 
	 * @return the {@link DifferentiableStatisticalModel} that is internally used
	 * 
	 * @throws CloneNotSupportedException
	 *             if the {@link DifferentiableStatisticalModel} could not be
	 *             cloned
	 */
	public DifferentiableStatisticalModel getFunction() throws CloneNotSupportedException {
		return (DifferentiableStatisticalModel)nsf.clone();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getNumberOfMotifs()
	 */
	public int getNumberOfMotifs() {
		if( nsf instanceof MotifDiscoverer ) {
			return ( (MotifDiscoverer)nsf ).getNumberOfMotifs();
		} else {
			return 0;
		}
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MutableMotifDiscoverer#adjustHiddenParameters(int, de.jstacs.data.DataSet[], double[][])
	 */
	public void adjustHiddenParameters( int index, DataSet[] data, double[][] weights ) throws Exception {
		if( nsf instanceof MutableMotifDiscoverer ) {
			DataSet[] mappedData = getMappedData( data );
			((MutableMotifDiscoverer)nsf).adjustHiddenParameters( index, mappedData, weights );
			for( int i = 0; i < data.length; i++ ) {
				mappedData[i] = null;
			}
			mappedData = null;
			System.gc();
		}
		
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MutableMotifDiscoverer#initializeMotif(int, de.jstacs.data.DataSet, double[])
	 */
	public void initializeMotif( int motifIndex, DataSet data, double[] weights ) throws Exception {
		if( nsf instanceof MutableMotifDiscoverer ) {
			DataSet[] mappedData = getMappedData( data );
			((MutableMotifDiscoverer)nsf).initializeMotif( motifIndex, mappedData[0], weights );
			mappedData[0] = null;
			mappedData = null;
			System.gc();
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MutableMotifDiscoverer#initializeMotifRandomly(int)
	 */
	public void initializeMotifRandomly( int motif ) throws Exception {
		if( nsf instanceof MutableMotifDiscoverer ) {
			((MutableMotifDiscoverer)nsf).initializeMotifRandomly( motif );
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MutableMotifDiscoverer#modifyMotif(int, int, int)
	 */
	public boolean modifyMotif( int motifIndex, int offsetLeft, int offsetRight ) throws Exception {
		if( nsf instanceof MutableMotifDiscoverer ) {
			return ((MutableMotifDiscoverer)nsf).modifyMotif( motifIndex, offsetLeft, offsetRight );
		} else {
			return false;
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getGlobalIndexOfMotifInComponent(int, int)
	 */
	public int getGlobalIndexOfMotifInComponent( int component, int motif ) {
		if( nsf instanceof MotifDiscoverer ) {
			return ((MotifDiscoverer)nsf).getGlobalIndexOfMotifInComponent( component, motif );
		} else {
			throw new ArrayIndexOutOfBoundsException();
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getIndexOfMaximalComponentFor(de.jstacs.data.Sequence)
	 */
	public int getIndexOfMaximalComponentFor( Sequence sequence ) throws Exception {
		if( nsf instanceof MotifDiscoverer ) {
			mappedSeq.setSequence( sequence );
			return ((MotifDiscoverer)nsf).getIndexOfMaximalComponentFor( mappedSeq );
		} else {
			throw new IllegalArgumentException();
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getMotifLength(int)
	 */
	public int getMotifLength( int motif ) {
		if( nsf instanceof MotifDiscoverer ) {
			return ((MotifDiscoverer)nsf).getMotifLength( motif );
		} else {
			throw new ArrayIndexOutOfBoundsException();
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getNumberOfComponents()
	 */
	public int getNumberOfComponents() {
		if( nsf instanceof MotifDiscoverer ) {
			return ((MotifDiscoverer)nsf).getNumberOfComponents();
		} else {
			return 1;
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getNumberOfMotifsInComponent(int)
	 */
	public int getNumberOfMotifsInComponent( int component ) {
		if( nsf instanceof MotifDiscoverer ) {
			return ((MotifDiscoverer)nsf).getNumberOfMotifsInComponent( component );
		} else {
			if( component > 0 ) {
				throw new IllegalArgumentException();
			} else {
				return 1;
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getProfileOfScoresFor(int, int, de.jstacs.data.Sequence, int, de.jstacs.motifDiscovery.MotifDiscoverer.KindOfProfile)
	 */
	public double[] getProfileOfScoresFor( int component, int motif, Sequence sequence, int startpos, KindOfProfile kind ) throws Exception {
		if( nsf instanceof MotifDiscoverer ) {
			mappedSeq.setSequence( sequence );
			double[] res = ((MotifDiscoverer)nsf).getProfileOfScoresFor( component, motif, mappedSeq, startpos, kind );
			double norm = mappedSeq.getLogNumberOfPossibleOriginalSequences( startpos, length == 0 ? sequence.getLength() : ( startpos + length ) );
			for( int i = 0; i < res.length; i++ ) {
				res[i] -= norm;
			}
			return res;
		} else {
			throw new IllegalArgumentException();
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getStrandProbabilitiesFor(int, int, de.jstacs.data.Sequence, int)
	 */
	public double[] getStrandProbabilitiesFor( int component, int motif, Sequence sequence, int startpos ) throws Exception {
		if( nsf instanceof MotifDiscoverer ) {
			mappedSeq.setSequence( sequence );
			return ((MotifDiscoverer)nsf).getStrandProbabilitiesFor( component, motif, mappedSeq, startpos );
		} else {
			throw new IllegalArgumentException();
		}
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.Mutable#modify(int, int)
	 */
	public boolean modify( int offsetLeft, int offsetRight ) {
		if( nsf instanceof Mutable ) {
			return ((Mutable) nsf).modify( offsetLeft, offsetRight );
		} else {
			return false;
		}
	}

	/**
	 * This class allows to change the original sequence. Be careful with this
	 * option since it might cause some trouble.
	 * 
	 * @author Jens Keilwagen
	 */
	private static final class MutableMappedDiscreteSequence extends MappedDiscreteSequence {

		private MutableMappedDiscreteSequence( AlphabetContainer originalAlphabetContainer, DiscreteAlphabetMapping[] transformation )
																																				throws WrongAlphabetException {
			super( originalAlphabetContainer, null, transformation );
		}

		private void setSequence( Sequence original ) throws IllegalArgumentException {
			if( original == null || !(original instanceof SimpleDiscreteSequence) || !originalAlphabetContainer.checkConsistency( original.getAlphabetContainer() ) ) {
				throw new IllegalArgumentException();
			}
			this.original = (SimpleDiscreteSequence) original;
		}

		private DiscreteAlphabetMapping[] getTransformations() {
			return transformation.clone();
		}

		private SimpleDiscreteSequence getOriginalSequence() {
			return original;
		}
	}
}
