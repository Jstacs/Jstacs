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

package de.jstacs.parameters;


/**
 * Interface for {@link Parameter}s that can be converted to and extracted from
 * <a href="http://galaxy.psu.edu/">Galaxy</a> representations.
 * The methods of this interface are used in the {@link de.jstacs.tools.ui.galaxy.GalaxyAdaptor} to create
 * Galaxy representation, i.e., config-files, from {@link ParameterSet}s that contain
 * the parameters for a specific application.
 * @author Jan Grau
 *
 */
public interface GalaxyConvertible {

	/**
	 * Creates an Galaxy XML-representation of the parameters and appends it to <code>descBuffer</code>
	 * and variable configuration and appends it to <code>configBuffer</code>. The variable configuration
	 * is also used to parse user-supplied values returned by Galaxy.
	 * 
	 * @param namePrefix the prefix of the variable name used in Galaxy
	 * @param configPrefix the prefix for conditionals
	 * @param depth the depth in the parameter hierarchy, used for graphical representation of nesting
	 * @param descBuffer the buffer for the parameter description
	 * @param configBuffer the buffer for the configuration line
	 * @param addLine if true, a line is added before the title of a parameter
	 * @param indentation the number of tabs that is used for indentation, if smaller than zero no indentation is used at all
	 * @throws Exception if the conversion fails
	 */
	public void toGalaxy(String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer, boolean addLine, int indentation ) throws Exception;
	
	/**
	 * Parses the contents of <code>command</code> in the format defined by <code>configBuffer</code> of {@link GalaxyConvertible#toGalaxy(String, String, int, StringBuffer, StringBuffer, boolean, int)}
	 * and sets the values of the {@link Parameter} or {@link ParameterSet} accordingly.
	 * @param namePrefix the prefix of the variable name
	 * @param command the command string
	 * @throws Exception if the command string could not be parsed
	 */
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception;
	
}
