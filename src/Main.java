import java.io.*;
import java.util.Scanner;
import java.util.ArrayList;

public class Main {

    // 1.1
    //je m'en criss un peu ici so j'ai tout criss dans un tableau, dans le fond
    //les index 0,2,4,6... sont les >dsad|dasd|sdads| et les index 1,3,5,... sont
    //les nucleotides, good luck
    //trop de trouble de les filtrer
    public static ArrayList<String> readFasta(String x) throws FileNotFoundException {

        ArrayList<String> database = new ArrayList<String>();
        try (Scanner sc = new Scanner(new File(x))) {
            while (sc.hasNextLine()) {
                database.add(sc.nextLine().trim());
            }
        }
        return database;
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

    //1.3
    //1.4
    //1.5
    //1.6
    //1.7
    public static void main(String[] args) throws FileNotFoundException {
        System.out.println(readFasta("src/tRNAs.fasta"));
        ArrayList<String> database = readFasta("src/tRNAs.fasta");
        ArrayList<String> unknown = readFasta("src/unknown.fasta");
        System.out.println(kmer(unknown.get(1), 11));
    }
}
