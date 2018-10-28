package projects.xanthogenomes;

import java.text.NumberFormat;

import de.jstacs.algorithms.alignment.Alignment;
import de.jstacs.algorithms.alignment.StringAlignment;
import de.jstacs.algorithms.alignment.Alignment.AlignmentType;
import de.jstacs.algorithms.alignment.PairwiseStringAlignment;
import de.jstacs.algorithms.alignment.cost.Costs;


public class TALEAligner {

	
	public static StringAlignment align(TALE tale1, TALE tale2, Costs costs, AlignmentType at, double extraGapOpening, double extraGapExtension){
		
		TALESequence s1 = new TALESequence( tale1 );
		TALESequence s2 = new TALESequence( tale2 );
		
		Alignment al = new Alignment( costs );
		
		String str1 = null;
		
		StringAlignment psa = null;
		if(s1.getLength() >= s2.getLength()){
			psa = al.getAlignment( at , s1, s2 );
			str1 = psa.getAlignedString( 1 );
		}else if(s1.getLength() < s2.getLength()){
			psa = al.getAlignment( at , s2, s1 );
			str1 = psa.getAlignedString( 1 );
			psa = new StringAlignment( psa.getCost(), psa.getAlignedString( 1 ), psa.getAlignedString( 0 ) );
			
		}
		
		
		if(at == AlignmentType.SEMI_GLOBAL){
			str1 = str1.replaceAll( " ", "" );
			double cost = psa.getCost();			
			
			if(str1.startsWith( "-" )){
				cost += extraGapOpening;
				for(int i=0;str1.charAt( i ) == '-';i+=2){
					cost += extraGapExtension;
				}
			}
			if(str1.endsWith( "-" )){
				cost += extraGapOpening;
				for(int i=str1.length()-1;str1.charAt( i ) == '-';i-=2){
					cost += extraGapExtension;
				}
			}
			psa = new StringAlignment( cost, psa.getAlignedString( 0 ), psa.getAlignedString( 1 ) );
		}
		
		return psa;
	}	
	
	
	public static String alignmentToString( StringAlignment sa, NumberFormat nf ) {
		String s1 = sa.getAlignedString( 0 );
		String s2 = sa.getAlignedString( 1 );
		StringBuffer pat = new StringBuffer();
		for(int i=0;i<s1.length();i++){
			if(s1.charAt( i ) == ' ' || s1.charAt( i ) == '-' || s2.charAt( i ) == '-'){
				pat.append( " " );
			}else{
				if(s1.charAt( i ) == s2.charAt( i )){
					pat.append( '|' );
				}else{
					pat.append( ':' );
				}
			}
		}
		return s1+"\n"+pat.toString()+"\n"+s2+"\nCost: "+nf.format( sa.getCost() );
		
	}
	
}
