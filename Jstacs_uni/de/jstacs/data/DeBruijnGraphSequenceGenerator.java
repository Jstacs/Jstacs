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

package de.jstacs.data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Random;

import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.CyclicSequenceAdaptor;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;

/**
 * Class for creating De Bruin sequences using explicit De Bruijn graphs.
 * 
 * @author Jan Grau
 *
 */
public class DeBruijnGraphSequenceGenerator {

	private static class Node{

		private String label;
		private ArrayList<Node> children;
		private ArrayList<Integer> edgeVisited;

		private int nIncoming;

		public Node(String label){
			this.label = label;
			this.children = new ArrayList<Node>();
			this.edgeVisited = new ArrayList<Integer>();
			this.nIncoming = 0;
		}

		public String toString(){
			return label;
		}

		public boolean addChild(Node end) {
			if(this.children.contains(end)){

				return false;
			}else{
				children.add(end);

				edgeVisited.add(0);
				end.nIncoming++;
				return true;
			}
		}

		public int getNOutgoing() {
			return children.size();
		}

		public int getNIncoming() {
			return nIncoming;
		}

		public String getLabel() {
			return label;
		}

		public LinkedList<Node> getUnvisitedEdges(){
			LinkedList<Node> edges = new LinkedList<Node>();
			for(int i=0;i<children.size();i++){
				if(edgeVisited.get(i) == 0){
					edges.add(children.get(i));
				}
			}
			return edges;
		}

		public Node getUnvisitedEdge(Random r) {
			/*for(int i=0;i<children.size();i++){
				if(edgeVisited.get(i) == 0){
					return children.get(i);
				}
			}
			return null;*/
			LinkedList<Node> edges = getUnvisitedEdges();
			if(edges.size() == 0){
				return null;
			}else{
				return edges.get(r.nextInt(edges.size()));
			}
		}

		public void setVisited(Node node){
			int i = children.indexOf(node);
			edgeVisited.set(i, edgeVisited.get(i)+1);
		}

	}



	/**
	 * Generates a De Bruijn sequence of length {@latex.inline $|A|^n$}, where A denotes the alphabet.
	 * @param alphabet the alphabet
	 * @param n the exponent of length length, corresponds to the length of n-mers covered exactly once
	 * @return the sequence (wrapped in an array)
	 * @throws WrongAlphabetException if the alphabet is 
	 * @throws WrongSequenceTypeException
	 */
	public static CyclicSequenceAdaptor[] generate(DiscreteAlphabet alphabet, int n) throws WrongAlphabetException, WrongSequenceTypeException {

		return new CyclicSequenceAdaptor[]{generate(alphabet, n, 0)};
	}

	/**
	 * Generates a De Bruijn sequence using the supplied alphabet and the given alphabet shift, i.e., for a cyclic shift of the symbols 
	 * of the alphabet.
	 * @param alphabet the alphabet
	 * @param n the length of the covered n-mers
	 * @param alphabetShift the alphabet shift (0 equals no shift)
	 * @return the De Bruijn sequence
	 * @throws WrongAlphabetException 
	 * @throws IllegalArgumentException 
	 */
	public static CyclicSequenceAdaptor generate(DiscreteAlphabet alphabet, int n, int alphabetShift) throws IllegalArgumentException, WrongAlphabetException {

		Random r = new Random(117);
		
		DiscreteSequenceEnumerator en = new DiscreteSequenceEnumerator(new AlphabetContainer(alphabet), n, false);

		HashMap<String, Node> nodes = new HashMap<String, Node>();

		String first = null;
		while(en.hasMoreElements()){
			String seq = en.nextElement().toString();
			
			String s1 = seq.substring(0, seq.length()-1);
			if(first == null){
				first = s1;
			}
			String s2 = seq.substring(1);

			if(!nodes.containsKey(s1)){
				nodes.put(s1, new Node(s1));
			}
			if(!nodes.containsKey(s2)){
				nodes.put(s2, new Node(s2));
			}

			Node n1 = nodes.get(s1);
			Node n2 = nodes.get(s2);
			n1.addChild(n2);
		}

		String str = findEulerPath(r, nodes.get(first));
		
		return new CyclicSequenceAdaptor(Sequence.create(new AlphabetContainer(alphabet), str));
		
	}

	private static String findEulerPath(Random r, Node start){


		LinkedList<Node> stack = new LinkedList<Node>();
		stack.addFirst(start);

		StringBuffer sb = new StringBuffer();

		while(stack.size() > 0){
			//System.out.println("before: "+stack);
			Node curr = stack.peek();
			Node next = null;
			if( (next = curr.getUnvisitedEdge(r) ) != null ){
				curr.setVisited(next);
				stack.addFirst(next);
				//System.out.println("add: "+next);
			}else{
				Node removed = stack.pop();
				
				sb.append(removed.getLabel().charAt(0));
			}
		}

		/*String str = sb.reverse().toString();
		String[] parts = str.split("\\$");

		sb = new StringBuffer();
		for(int i=parts.length-1;i>=0;i--){
			if(i==0){
				parts[i] = parts[i].substring(1);
			}
			sb.append(parts[i]);
		}*/

		return sb.toString().substring(0, sb.length()-1);
	}



}
