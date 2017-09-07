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

package de.jstacs.results.savers;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import de.jstacs.results.Result;

/**
 * Library of all registered {@link ResultSaver}s.
 * 
 * @author Jan Grau
 *
 */
public class ResultSaverLibrary {

	//registers standard ResultSavers
	static{
		TextResultSaver.register();
		ResultSetResultSaver.register();
		ListResultSaver.register();
		DataSetResultSaver.register();
		PlotGeneratorResultSaver.register();
		StorableResultSaver.register();
	}
	
	private static HashMap<Class<? extends Result>, ResultSaver> map;
	
	/**
	 * Registers the supplied {@link ResultSaver} for the given class in the library.
	 * @param clazz the class
	 * @param saver the {@link ResultSaver}
	 * @param <T> the class of the result
	 */
	public static <T extends Result> void register(Class<? extends T> clazz, ResultSaver<T> saver){
		if(map == null){
			map = new HashMap<Class<? extends Result>, ResultSaver>();
		}
		map.put( clazz, saver );
	}
	
	/**
	 * Returns the most suitable {@link ResultSaver} (if any) currently registered in the library.
	 * If a {@link ResultSaver} for the given {@link Result} type is registered, this one is returned.
	 * Otherwise, the list of {@link ResultSaver}s is searched for a {@link ResultSaver} registered for
	 * a superclass of <code>T</code> and the first hit is returned.
	 * 
	 * @param resultClass the class of the result to be saved
	 * @param <T> the class of the result
	 * @return the appropriate {@link ResultSaver} or <code>null</code> if not such {@link ResultSaver} has been registered
	 */
	public static <T extends Result> ResultSaver<T> getSaver(Class<? extends T> resultClass){
		
		ResultSaver ren = map.get( resultClass );
		if(ren == null){
			Set<Class<? extends Result>> clazzes = map.keySet();
			Iterator<Class<? extends Result>> it = clazzes.iterator();
			while(it.hasNext()){
				Class<? extends Result> clazz = it.next();
				if(clazz.isAssignableFrom( resultClass )){
					return map.get( clazz );
				}
			}
		}
		if(ren == null){
			System.err.println( "Did not find a saver for "+resultClass+". Custom savers need to be registered by "+ResultSaverLibrary.class.getName()+".register()." );
		}
		return ren;
	}

}
