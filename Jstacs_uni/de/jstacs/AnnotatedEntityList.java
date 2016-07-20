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

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;

/**
 * Class for a list of {@link AnnotatedEntity}s where
 * elements can be accessed either by index or by the name 
 * of the {@link AnnotatedEntity}. 
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
 * @param <T> an extension of {@link AnnotatedEntity} 
 */
public class AnnotatedEntityList<T extends AnnotatedEntity> {
	
	private ArrayList<T> entityList;
	private HashMap<String, T> entities;
	
	/**
	 * Creates a new {@link AnnotatedEntityList} with an initial
	 * capacity of 10.
	 */
	public AnnotatedEntityList() {
		this( 10 );
	}
	
	/**
	 * Creates a new {@link AnnotatedEntityList} with given initial
	 * capacity.
	 * 
	 * @param initialCapacity the initial capacity
	 */
	public AnnotatedEntityList( int initialCapacity ) {
		entityList = new ArrayList<T>(initialCapacity);
		entities = new HashMap<String, T>(initialCapacity);
	}
	
	/**
	 * Replaces the {@link AnnotatedEntity} at index <code>idx</code> with
	 * the {@link AnnotatedEntity} <code>entity</code>
	 * @param idx the index
	 * @param entity the new {@link AnnotatedEntity} at index <code>idx</code>
	 */
	public void set( int idx, T entity ){
		String name = entityList.get( idx ).getName();
		entities.remove( name );
		entityList.set( idx, entity );
		entities.put( entity.getName(), entity );
	}
	
	/**
	 * Removes and returns the {@link AnnotatedEntity} at index <code>idx</code>.
	 * @param idx the index of the {@link AnnotatedEntity} to be removed
	 * @return the {@link AnnotatedEntity} that has been at index <code>idx</code>
	 */
	public T remove(int idx){
		T p = entityList.remove( idx );
		entities.remove( p.getName() );
		return p;
	}
	
	/**
	 * Adds the {@link AnnotatedEntity} <code>entity</code>
	 * at index <code>idx</code> to the list. If <code>entities</code> contains {@link AnnotatedEntity}s with duplicate names,
	 * or an {@link AnnotatedEntity} with identical name already exists in the list,
	 * an {@link IllegalArgumentException} is thrown.
	 * @param idx the index
	 * @param entity the added {@link AnnotatedEntity}
	 * @see java.util.List#add(int, Object)
	 */
	public void add( int idx, T entity ){
		if( entities.containsKey( entity.getName() ) ) {
			throw new IllegalArgumentException( "The name \"" +entity.getName()+ "\" is already contained." );
		}
		entityList.add( idx, entity );
		entities.put( entity.getName(), entity );
	}
	
	/**
	 * Adds all {@link AnnotatedEntity}s in <code>entities</code>
	 * to the list. If <code>entities</code> contains {@link AnnotatedEntity}s with duplicate names,
	 * or an {@link AnnotatedEntity} with identical name already exists in the list,
	 * an {@link IllegalArgumentException} is thrown.
	 * @param entities the added {@link AnnotatedEntity}
	 */
	public void add( T... entities ) {
		for(int i=0;i<entities.length;i++){
			if( this.entities.containsKey( entities[i].getName() ) ) {
				throw new IllegalArgumentException( "The name \"" +entities[i].getName()+ "\" is already contained." );
			}
			entityList.add( entities[i] );
			this.entities.put( entities[i].getName(), entities[i] );
		}
	}
	
	/**
	 * Adds all {@link AnnotatedEntity}s in <code>entities</code>
	 * to the list. If <code>entities</code> contains {@link AnnotatedEntity}s with duplicate names,
	 * or an {@link AnnotatedEntity} with identical name already exists in the list,
	 * an {@link IllegalArgumentException} is thrown.
	 * 
	 * @param entities the added {@link AnnotatedEntity}
	 */
	@SuppressWarnings("unchecked")
	public void addAll(Collection<? extends T> entities){
		Iterator<? extends T> it = entities.iterator();
		while(it.hasNext()){
			this.add( it.next() );
		}
	}
	
	/**
	 * Returns the {@link AnnotatedEntity} at index <code>index</code>
	 * in the list.
	 * @param index the index
	 * @return the {@link AnnotatedEntity} at <code>index</code>
	 */
	public T get( int index ) {
		return entityList.get(index);
	}

	/**
	 * Returns the {@link AnnotatedEntity} with name <code>name</code>
	 * in the list.
	 * @param name the name
	 * @return the {@link AnnotatedEntity} with name <code>name</code>
	 * @see AnnotatedEntity#getName()
	 */
	public T get( String name ) {
		return entities.get( name );
	}
	
	/**
	 * Returns the number of {@link AnnotatedEntity}s (not the capacity)
	 * in the {@link AnnotatedEntityList}.
	 * @return the number of {@link AnnotatedEntity}s
	 */
	public int size() {
		return entityList.size();
	}
	
	/**
	 * Returns the names of all {@link AnnotatedEntity}s in the list.
	 * @return the names
	 * @see AnnotatedEntity#getName()
	 */
	public String[] getNames(){
		String[] names = new String[entityList.size()];
		for( int i = 0; i < names.length; i++ ) {
			names[i] = entityList.get(i).getName();
		}
		return names; //sorted
		//return entities.keySet().toArray( new String[0] );//not sorted
	}
	
	/**
	 * Returns the {@link AnnotatedEntity}s in this list
	 * as an array. This method behaves exactly like {@link ArrayList#toArray(Object[])}.
	 * @param <E> a sub-type of {@link AnnotatedEntity}
	 * @param ar an array of the desired type
	 * @return the {@link AnnotatedEntity}s in this list as an array
	 */
	public <E extends T> E[] toArray(E[] ar){
		return entityList.toArray(ar);
	}
}
