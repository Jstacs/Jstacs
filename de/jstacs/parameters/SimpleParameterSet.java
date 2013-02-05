/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.parameters;

import de.jstacs.io.NonParsableException;

/**
 * Class for a {@link ParameterSet} that is constructed from an array of {@link Parameter}s.
 * 
 * @author Jan Grau
 */
public class SimpleParameterSet extends ParameterSet{

	/**
	 * Creates a new <code>SimpleParameterSet</code> from an array of <code>Parameter</code>s.
	 * @param parameters the parameters
	 */
	public SimpleParameterSet(Parameter... parameters){
		super(parameters);
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link SimpleParameterSet} from its XML representation.
	 * 
	 * @param representation the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException if the {@link StringBuffer} <code>representation</code> could not be parsed
	 */
	public SimpleParameterSet(StringBuffer representation) throws NonParsableException{
		super(representation);
	}
	
	/* (non-Javadoc)
	 * @see de.jstacs.parameters.ParameterSet#clone()
	 */
	@Override
	public SimpleParameterSet clone() throws CloneNotSupportedException{
		SimpleParameterSet clone = (SimpleParameterSet)super.clone();
		return clone;
	}
    
    /* (non-Javadoc)
     * @see de.jstacs.parameters.ParameterSet#reset()
     */
    @Override
    public void reset(){
        if( parameters == null )
        {
			return;
        }
        for( int i = 0; i < parameters.size(); i++ )
        {
            parameters.get( i ).reset();
        }
    }
    
    
    public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer, boolean[] addLine ) throws Exception {
		for(int i=0;i<getNumberOfParameters();i++){
			((GalaxyConvertible)getParameterAt( i )).toGalaxy( namePrefix+"_ps", configPrefix, depth+1, descBuffer, configBuffer, addLine == null ? false : addLine[i] );
			descBuffer.append( "\n" );
			configBuffer.append( "\n" );
		}	
	}
}
