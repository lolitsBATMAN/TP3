import java.io.*;
import java.util.Scanner;
import java.util.ArrayList;

public class Main {
    private ArrayList<DBentry> database = new ArrayList<DBentry>();

    // 1.1
    //je m'en criss un peu ici so j'ai tout criss dans un tableau, dans le fond
    //les index 0,2,4,6... sont les >dsad|dasd|sdads| et les index 1,3,5,... sont
    //les nucleotides, good luck
    //trop de trouble de les filtrer

    //2e argument la dabatase quon creer dans la main
    public void readFasta(String file) throws FileNotFoundException {


        try (Scanner sc = new Scanner(new File(file))) {

            while (sc.hasNextLine()) {
                String infos = sc.nextLine();
                String sequence = sc.nextLine();

                //create entry that contains these 2 strings
                DBentry currEntry = new DBentry(infos,sequence);

                //add entry to database
                this.database.add(currEntry);
            }
        }
    }

    //1.2
    //la position de chaque KMER est l'index dans le tableau
    public static ArrayList<String> kmer (String word, int k){
        ArrayList<String> mer = new ArrayList<String>();

        for (int i  = 0; i <= word.length() - k ; i++){
            mer.add(word.substring(i, i+k));
        }

        return mer;
    }

    public static ArrayList<String> hsp(String input, String data, String seed){
        int k = seed.length();
        ArrayList<String> hsps = new ArrayList<String>();
        ArrayList<String> kinput = kmer(input, k);
        ArrayList<String> kdata = kmer(data, k);

        for(int i = 0; i < kinput.size(); i++){
            for(int j = 0; j < kdata.size(); j++){
                Boolean bruh = true;
                for (int z = 0; z < k; z++){
                    if (kinput.get(i).charAt(z) != kdata.get(j).charAt(z) && seed.charAt(z) == '1') {
                        bruh = false;
                        break;
                    } else if(kinput.get(i).charAt(z) == kdata.get(j).charAt(z) && seed.charAt(z) == '0') {
                        bruh = false;
                        break;
                    }
                }
                if (bruh){
                    hsps.add(kinput.get(i));
                }
            }
        }

        return hsps;
    }

    //1.3
    //1.4
    //1.5
    //1.6
    //1.7
    public static void main(String[] args) throws FileNotFoundException {
        //System.out.println(readFasta("src/tRNAs.fasta"));
        //ArrayList<DBentry> database = readFasta("src/tRNAs.fasta");
        //ArrayList<String> unknown = readFasta("src/unknown.fasta");
        //System.out.println(kmer(unknown.get(1), 11));

        //Creer ici la database pour avoir acces partout
        System.out.println(hsp("BAAA","AAAB" , "111"));
    }
}
