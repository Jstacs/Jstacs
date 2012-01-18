package de.jstacs.classifiers.neuralNetwork.neurons;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;


public class InputNeuron extends Neuron {

	private int position;
	
	public InputNeuron(int position){
		super();
		this.position = position;
	}
	
	@Override
	public int getNumberOfWeights() {
		return 0;
	}

	@Override
	public double getOutput( Sequence input ) {
		return input.continuousVal( position );
	}

	@Override
	public void initializeRandomly() {		
	}

	@Override
	public void reset() {
	}

	@Override
	public double getError( Sequence input, double weight, double[] desiredOutputs ) {
		double err = 0;
		for(int i=0;i<getNumberOfDescendants();i++){
			err += getDescendant( i ).getError( input, weight, desiredOutputs );
		}
		return err;
	}

	@Override
	public void adaptWeights( double eta ) {
	}



	@Override
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
		position = XMLParser.extractObjectForTags( xml, "position", int.class );
	}



	@Override
	protected StringBuffer getFurtherInformation( ) {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, position, "position" );
		return xml;
	}

	

}
