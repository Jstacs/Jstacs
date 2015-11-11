package de.jstacs.data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Random;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.CyclicSequenceAdaptor;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;


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



	public static CyclicSequenceAdaptor[] generate(DiscreteAlphabet alphabet, int n) throws WrongAlphabetException, WrongSequenceTypeException, OperationNotSupportedException{

		return new CyclicSequenceAdaptor[]{generate(alphabet, n, 0)};
	}

	public static CyclicSequenceAdaptor generate(DiscreteAlphabet alphabet, int n, int alphabetShift) throws WrongAlphabetException, WrongSequenceTypeException{

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
