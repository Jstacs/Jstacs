package de.jstacs.data.sequences;

import javax.naming.OperationNotSupportedException;

/**
 * This class is an adapter for sequence to mimic cyclic sequences.
 * 
 * @author Jens Keilwagen, Jan Grau
 *
 * @param <T> the type of the internal sequence
 */
public class CyclicSequenceAdaptor<T> extends Sequence<T> {

	//TODO wrap return values of rc, comp
	
	private Sequence<T> seq;
	private int l;
	private int extLength;
	
	/**
	 * Creates a new cyclic sequence of given virtual length (i.e., the length reported by {@link #getLength()}).
	 * @param seq the original sequence
	 * @param extLength the virtual length
	 */
	public CyclicSequenceAdaptor(Sequence<T> seq, int extLength) {
		super(seq.getAlphabetContainer(), seq.getAnnotation());
		this.seq = seq;
		this.l = this.seq.getLength();
		this.extLength = extLength;
	}
	
	/**
	 * Creates a new cyclic sequence of the length of the original sequence.
	 * @param seq the original sequence
	 */
	public CyclicSequenceAdaptor(Sequence<T> seq) {
		this(seq,seq.getLength());
	}

	@Override
	public int compareTo(Sequence<T> s) {
		// TODO cyclic?
		return seq.compareTo(s);
	}

	@Override
	public double continuousVal(int pos) {
		return seq.continuousVal(pos%l);
	}

	@Override
	public int discreteVal(int pos) {
		return seq.discreteVal(pos%l);
	}

	@Override
	protected CyclicSequenceAdaptor<T> flatCloneWithoutAnnotation() {
		// TODO ??
		return new CyclicSequenceAdaptor<T>(seq);
	}

	@Override
	public int getLength() {
		return extLength;
	}

	@Override
	protected int compareTo(T t1, T t2) {
		return seq.compareTo(t1, t2);
	}

	@Override
	protected Object getEmptyRepresentation() {
		return seq.getEmptyRepresentation();
	}

	@Override
	protected void addToRepresentation(Object representation, int pos, String delim) {
		seq.addToRepresentation(representation, pos%l, delim);
	}

	@Override
	protected String getStringRepresentation(Object representation) {
		// TODO cyclic?
		return seq.getStringRepresentation(representation);
	}

	@Override
	protected int hashCodeForPos(int pos) {
		return seq.hashCodeForPos(pos%l);
	}

	@Override
	public boolean isMultiDimensional() {
		return seq.isMultiDimensional();
	}

	@Override
	public T getEmptyContainer() {
		return seq.getEmptyContainer();
	}

	@Override
	public void fillContainer(T container, int pos) {
		seq.fillContainer(container, pos%l);
	}

	@Override
	public CyclicSequenceAdaptor<T> reverse( int start, int end ) throws OperationNotSupportedException {
		return new CyclicSequenceAdaptor<T>( seq.reverse( start, end ) );
	}

	@Override
	public CyclicSequenceAdaptor<T> complement() throws OperationNotSupportedException {
		return new CyclicSequenceAdaptor<T>( seq.complement() );
	}

	@Override
	public CyclicSequenceAdaptor<T> reverseComplement() throws OperationNotSupportedException {
		return new CyclicSequenceAdaptor<T>( seq.reverseComplement() );
	}

	@Override
	public CyclicSequenceAdaptor<T> complement( int start, int end ) throws OperationNotSupportedException {
		return new CyclicSequenceAdaptor<T>( seq.complement( start, end ) );
	}

	@Override
	public CyclicSequenceAdaptor<T> reverseComplement( int start, int end ) throws OperationNotSupportedException {
		return new CyclicSequenceAdaptor<T>( seq.reverseComplement( start, end ) );
	}
	
	/**
	 * Returns a new cyclic sequence using the internal sequence of this {@link CyclicSequenceAdaptor} but with 
	 * the supplied virtual length
	 * @param length the virtual length
	 * @return the cyclic sequence
	 * @see #CyclicSequenceAdaptor(Sequence, int)
	 */
	public CyclicSequenceAdaptor<T> getSuperSequence(int length){
		return new CyclicSequenceAdaptor<T>( seq, length );
	}
	
	
	
}