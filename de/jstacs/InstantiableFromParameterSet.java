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

package de.jstacs;

import de.jstacs.parameters.InstanceParameterSet;

/**
 * Interface for all classes that can be instantiated from a
 * {@link InstanceParameterSet}. Each class <code>ABC</code> implementing this
 * interface <b>must</b> provide a constructor of kind<br>
 * <code>ABC({@link InstanceParameterSet} parameters)</code>.
 * 
 * <br>
 * <br>
 * 
 * In order to provide a suitable {@link InstanceParameterSet} to obtain all
 * necessary values of these variables the user should implement a subclass of
 * the abstract class {@link InstanceParameterSet} or one of its subclasses,
 * that creates the necessary {@link de.jstacs.parameters.Parameter}s in its
 * <code>loadParameters()</code>-method.
 * 
 * @author Jan Grau, Jens Keilwagen
 * @deprecated Parameters are used for user parameters of tools now but no longer for directly instantiating objects
 */
@Deprecated 
public interface InstantiableFromParameterSet {

	/**
	 * Returns the {@link InstanceParameterSet} that has been used to
	 * instantiate the current instance of the implementing class. If the
	 * current instance was not created using an {@link InstanceParameterSet},
	 * an equivalent {@link InstanceParameterSet} should be returned, so that an
	 * instance created using this {@link InstanceParameterSet} would be in
	 * principle equal to the current instance.
	 * 
	 * @return the current {@link InstanceParameterSet}
	 * 
	 * @throws Exception
	 *             if the {@link InstanceParameterSet} could not be returned
	 */
	public InstanceParameterSet<? extends InstantiableFromParameterSet> getCurrentParameterSet() throws Exception;

}