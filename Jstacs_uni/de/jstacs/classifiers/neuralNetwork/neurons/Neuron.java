package de.jstacs.classifiers.neuralNetwork.neurons;

import java.util.ArrayList;

import de.jstacs.Storable;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;


public abstract class Neuron implements Storable{
	
	private ArrayList<InnerNeuron> descendants;
	
	public Neuron(){
		this.descendants = new ArrayList<InnerNeuron>();
	}
	
	public Neuron(StringBuffer xml) throws NonParsableException{
		fromXML(xml);
		if(this.descendants == null){
			throw new NonParsableException();
		}
	}
	
	private void fromXML( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, getClass().getSimpleName() );
		extractFurtherInformation(xml);
	}

	protected abstract void extractFurtherInformation( StringBuffer xml ) throws NonParsableException;

	protected abstract StringBuffer getFurtherInformation( );

	public StringBuffer toXML(){
		StringBuffer sb = new StringBuffer();
		sb.append( getFurtherInformation( ));
		XMLParser.addTags( sb, getClass().getSimpleName() );
		return sb;
	}
	
	public final int getNumberOfDescendants(){
		return descendants.size();
	}
	
	public final InnerNeuron getDescendant(int i){
		return descendants.get( i );
	}
	
	protected final void addDescendant( InnerNeuron neuron ){
		this.descendants.add( neuron );
	}
		
	public abstract int getNumberOfWeights();
	
	public abstract double getOutput(Sequence input);
	
	public abstract void initializeRandomly();
	
	public abstract void reset();
	
	public abstract double getError(Sequence input, double weight, double[] desiredOutputs );

	public abstract void adaptWeights( double eta );
}
