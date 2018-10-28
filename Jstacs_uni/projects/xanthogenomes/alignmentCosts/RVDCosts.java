package projects.xanthogenomes.alignmentCosts;

import de.jstacs.algorithms.alignment.cost.Costs;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import projects.xanthogenomes.TALE;
import projects.xanthogenomes.TALESequence;


public class RVDCosts implements Costs {

	private double gap;
	private double twelve, thirteen, bonus;
	
	public RVDCosts(double gap, double twelve, double thirteen, double bonus){
		this.gap = gap;
		this.twelve = twelve;
		this.thirteen = thirteen;
		this.bonus = bonus;
	}
	
	public RVDCosts(StringBuffer xml) throws NonParsableException{
		xml = XMLParser.extractForTag( xml, "RVDCosts" );
		bonus = (Double)XMLParser.extractObjectForTags( xml, "bonus" );
		gap = (Double)XMLParser.extractObjectForTags( xml, "gap" );
		thirteen = (Double)XMLParser.extractObjectForTags( xml, "thirteen" );
		twelve = (Double)XMLParser.extractObjectForTags( xml, "twelve" );
	}
	
	@Override
	public double getCostFor( Sequence s1, Sequence s2, int i, int j ) {
		TALESequence t1 = (TALESequence)s1;
		TALESequence t2 = (TALESequence)s2;
		TALE tale1 = t1.getTALE();
		TALE tale2 = t2.getTALE();
		String rvd1 = tale1.getRepeat( i-1 ).getRvd();
		String rvd2 = tale2.getRepeat( j-1 ).getRvd();

		double cost = getCostFor( rvd1, rvd2 );
		
		return cost;
	}
	
	public double getCostFor(AlphabetContainer rvdAlp, int rvd1, int rvd2){
		return getCostFor( rvdAlp.getSymbol( 0, rvd1 ), rvdAlp.getSymbol( 0, rvd2 ) );
	}

	private double getCostFor(String rvd1, String rvd2){
		double cost = 0;
		if(rvd1.charAt(0) != rvd2.charAt(0)){
			cost += twelve;
		}
		if(rvd1.charAt(1) != rvd2.charAt(1)){
			cost += thirteen;
		}
		/*if(rvd1.charAt(0) != rvd2.charAt(0) || rvd1.charAt(1) != rvd2.charAt(1)){
			cost += twelve + thirteen;
		}*/
		if(cost==0){
			cost += bonus;
		}
		
		return cost;
	}
	
	
	public double getInsertCosts(){return getGapCosts();}
	
	public double getDeleteCosts(){return getGapCosts();}
	
	
	public double getGapCosts() {
		return gap;
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, bonus, "bonus" );
		XMLParser.appendObjectWithTags( xml, gap, "gap" );
		XMLParser.appendObjectWithTags( xml, thirteen, "thirteen" );
		XMLParser.appendObjectWithTags( xml, twelve, "twelve" );
		XMLParser.addTags( xml, "RVDCosts" );
		return xml;
	}

}
