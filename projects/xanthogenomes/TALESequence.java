package projects.xanthogenomes;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.sequences.Sequence;


public class TALESequence extends Sequence<TALE> {

	TALE internal;
	
	public TALESequence(TALE tale){
		super(tale.getRvdSequence().getAlphabetContainer(),null);
		this.internal = tale;
	}

	public TALE getTALE(){
		return internal;
	}
	
	@Override
	public double continuousVal( int pos ) {
		return internal.getRvdSequence().continuousVal( pos );
	}

	@Override
	public int discreteVal( int pos ) {
		return internal.getRvdSequence().discreteVal( pos );
	}

	@Override
	protected Sequence flatCloneWithoutAnnotation() {
		return null;
	}

	@Override
	public int getLength() {
		return internal.getNumberOfRepeats();
	}

	@Override
	protected int compareTo( TALE t1, TALE t2 ) {
		return 0;
	}

	@Override
	protected Object getEmptyRepresentation() {
		return null;
	}

	@Override
	protected void addToRepresentation( Object representation, int pos, String delim ) {
		
	}

	@Override
	protected String getStringRepresentation( Object representation ) {
		return "";
	}

	@Override
	protected int hashCodeForPos( int pos ) {
		return 0;
	}

	@Override
	public boolean isMultiDimensional() {
		return false;
	}

	@Override
	public TALE getEmptyContainer() {
		return null;
	}

	@Override
	public void fillContainer( TALE container, int pos ) {

	}

}
