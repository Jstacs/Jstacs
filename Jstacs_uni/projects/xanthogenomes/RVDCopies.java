package projects.xanthogenomes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Random;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.FileManager;
import projects.xanthogenomes.TALEFamilyBuilder.TALEFamily;

public class RVDCopies {

	private static class Match{
		
		private int start1, start2;
		private int len1, len2;
		private int lenm;
		
		public Match(int start1, int start2, Sequence rvds1, Sequence rvds2, int len){
			this.start1 = start1;
			this.start2 = start2;
			this.len1 = rvds1.getLength();
			this.len2 = rvds2.getLength();
			this.lenm = len;
		}
		
		
	}
	
	public static void main(String[] args) throws Exception {
		TALEFamilyBuilder builder = new TALEFamilyBuilder(FileManager.readFile(args[0]));
		
		TALEFamily[] fams = builder.getFamilies();

		int minlen = 4;
		
		
		
		
		Sequence[][] allSeqs = new Sequence[fams.length][];
		String[] names = new String[fams.length];
		String[] strains = new String[fams.length];
		
		for(int i=0;i<fams.length;i++){
			names[i] = fams[i].getFamilyId();
			allSeqs[i] = new Sequence[fams[i].getFamilySize()];
			String strain = null;
			for(int j=0;j<allSeqs[i].length;j++){
				allSeqs[i][j] = fams[i].getFamilyMembers()[j].getRvdSequence();
				strain = fams[i].getFamilyMembers()[j].getStrain();
				strain = strain.substring(0,strain.indexOf(" "));
				if(strains[i] == null){
					strains[i] = strain;
				}else if(!strains[i].equals(strain)){
					strains[i] = "mult";
				}
			}
		}
		
		
		AlphabetContainer rvdAlph = allSeqs[0][0].getAlphabetContainer(); 
		
		
		printCounts(allSeqs, names, "real",strains);
		
		System.exit(1);
		
		for(int r=0;r<1000;r++){
			Sequence[][] allSeqs2 = new Sequence[allSeqs.length][];
			for(int i=0;i<fams.length;i++){
				String[] al = fams[i].getInducedMultipleAlignment().getSecondElement();

				String[][] split = new String[al.length][];
				for(int k=0;k<al.length;k++){
					split[k] = al[k].trim().split(" ");
				}
				int[] part = new int[split[0].length];
				int curr = 0;
				for(int j=1;j<part.length;j++){
					innerloop:
						for(int k=0;k<split.length;k++){
							if( (split[k][j].equals("--") && !split[k][j-1].equals("--")) ||
									(!split[k][j].equals("--") && split[k][j-1].equals("--"))){
								curr++;
								break innerloop;
							}
						}
				part[j] = curr;
				}
				//	System.out.println(Arrays.toString(part));
				int[] ends = new int[part[part.length-1]+1];
				int j=-1,k=0;
				while(j<part.length){
					j++;
					while(j<part.length && (j==0 || part[j] == part[j-1]) ){
						j++;
					}
					//System.out.println(k+" "+j);
					ends[k] = j;
					k++;
				}
				//System.out.println(Arrays.toString(ends));
				int start = 0;
				int[] perms = new int[part.length];

				for(j=0;j<ends.length;j++){
					int end = ends[j];
					int[] perm = new int[end-start];
					for(int l=0;l<perm.length;l++){
						perm[l] = start+l;
					}
					permute(perm);

					System.arraycopy(perm, 0, perms, start, perm.length);

					start = ends[j];
				}
				//System.out.println(Arrays.toString(perms));


				if(split.length != allSeqs[i].length){
					throw new RuntimeException();
				}

				allSeqs2[i] = new Sequence[split.length];

				for(j=0;j<split.length;j++){
					String[] temp = split[j].clone();
					for(k=0;k<temp.length;k++){
						temp[k] = split[j][perms[k]];
					}
					StringBuffer strb = new StringBuffer();
					for(k=0;k<temp.length;k++){
						if(strb.length()>0){
							strb.append(" ");
						}
						if(!temp[k].equals("--")){
							strb.append(temp[k]);
						}
					}
					Sequence rvds = Sequence.create(rvdAlph, strb.toString(), " ");
					allSeqs2[i][j] = rvds;
					//System.out.println(allSeqs[i][j]+" <-> "+allSeqs2[i][j]);
				}

			}

			printCounts(allSeqs2, names, "rand"+r,strains);

		}
		
		
		
		
	}

	
	
	private static void printCounts(Sequence[][] allSeqs2, String[] names, String suffix,String[] strains){
		
		int minlen = 4;
		ArrayList<Match> matches = new ArrayList<RVDCopies.Match>();
		
		Sequence[][] allSeqs = new Sequence[allSeqs2.length][];
		
		for(int i=0;i<allSeqs2.length;i++){
			LinkedList<Sequence> seqs1 = new LinkedList<Sequence>();
			for(int j=0;j<allSeqs2[i].length;j++){
				Sequence rvds1 = allSeqs2[i][j];
				if(!seqs1.contains(rvds1)){
					seqs1.add(rvds1);
				}
			}
			allSeqs[i] = seqs1.toArray(new Sequence[0]);
		}
		
		
		
		
		for(int i=0;i<allSeqs.length;i++){
			Sequence[] seqs1 = allSeqs[i];
			
			for(int j=0;j<seqs1.length;j++){
				Sequence rvds1 = seqs1[j];
				for(int k=i+1;k<allSeqs.length;k++){
					Sequence[] seqs2 = allSeqs[k];
					for(int l=0;l<seqs2.length;l++){
						Sequence rvds2 = seqs2[l];
						matches.clear();
						findMatches(rvds1,rvds2,matches, minlen);
						
						for(int m=0;m<matches.size();m++){
							Match ma = matches.get(m);
							System.out.println(names[i]+"\t"+names[k]+"\t"
									+ ma.start1+"\t"+ma.start2+"\t"+(ma.len1-(ma.start1+ma.lenm))+"\t"+(ma.len2-(ma.start2+ma.lenm))+"\t"+ma.lenm+"\t"+suffix+"\t"+strains[i]+"\t"+strains[k]);
						}
						
					}
				}
				
			}
		}
		
	}
	
	
	private static void permute(int[] perm) {
		Random r = new Random();
		for(int i=0;i<perm.length;i++){
			int swap = i+r.nextInt(perm.length-i);
			int temp = perm[swap];
			perm[swap] = perm[i];
			perm[i] = temp;
		}
		
	}

	private static void findMatches(Sequence rvds1, Sequence rvds2, ArrayList<Match> matches, int minlen) {
		
		for(int i=0;i<rvds1.getLength()-minlen+1;i++){
			Sequence sub1 = rvds1.getSubSequence(i);
			for(int j=0;j<rvds2.getLength()-minlen+1;j++){
				Sequence sub2 = rvds2.getSubSequence(j);
				
				if(i==0 || j==0 || rvds1.discreteVal(i-1) != rvds2.discreteVal(j-1)){
					int k=0;
					while(k< sub1.getLength() && k < sub2.getLength() && sub1.discreteVal(k)==sub2.discreteVal(k)){
						k++;
					}
					if(k > minlen){
						matches.add(new Match(i,j,rvds1,rvds2,k));
					}
				}
				
				
			}
		}
		
	}

}
