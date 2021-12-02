import java.io.*;
import java.util.Scanner;
import java.util.ArrayList;

public class Main {
    // 1.1
    //je m'en criss un peu ici so j'ai tout criss dans un tableau, dans le fond
    //les index 0,2,4,6... sont les >dsad|dasd|sdads| et les index 1,3,5,... sont
    //les nucleotides, good luck
    //trop de trouble de les filtrer

    //2e argument la dabatase quon creer dans la main
    public static ArrayList<DBentry> readFasta(String file) throws FileNotFoundException {
        ArrayList<DBentry> database =  new ArrayList<DBentry>();
        try (Scanner sc = new Scanner(new File(file))) {

            while (sc.hasNextLine()) {
                String infos = sc.nextLine();
                String sequence = sc.nextLine();

                //create entry that contains these 2 strings
                DBentry currEntry = new DBentry(infos,sequence);

                //add entry to database
                database.add(currEntry);
            }
        }
        return database;
    }

    //1.2 and 1.3
    //la position de chaque KMER est l'index dans le tableau
    public static ArrayList<String> kmer (String word, int k){
        ArrayList<String> mer = new ArrayList<String>();

        for (int i  = 0; i <= word.length() - k ; i++){
            mer.add(word.substring(i, i+k));
        }

        return mer;
    }

    //retourne une list avec les hsps
    public static ArrayList<ArrayList<String>> hsp(String input, ArrayList<DBentry> data, String seed, ArrayList<ArrayList<Integer>> position, ArrayList<ArrayList<Integer>> positioninput){
        int k = seed.length();
        ArrayList<ArrayList<String>> hsps = new ArrayList<ArrayList<String>>();
        ArrayList<String> kinput = kmer(input, k);

        for (int y = 0; y < data.size(); y++) {
            ArrayList<String> kdata = kmer(data.get(y).getSequence(), k);
            hsps.add(new ArrayList<String>());
            position.add(new ArrayList<Integer>());
            positioninput.add(new ArrayList<Integer>());
            for (int i = 0; i < kinput.size(); i++) {
                for (int j = 0; j < kdata.size(); j++) {
                    Boolean bruh = true;
                    for (int z = 0; z < k; z++) {
                        if (kinput.get(i).charAt(z) != kdata.get(j).charAt(z) && seed.charAt(z) == '1') {
                            bruh = false;
                            break;
                        } else if (kinput.get(i).charAt(z) == kdata.get(j).charAt(z) && seed.charAt(z) == '0') {
                            bruh = false;
                            break;
                        }
                    }
                    if (bruh) {
                        hsps.get(y).add(kinput.get(i));
                        position.get(y).add(j);
                        positioninput.get(y).add(i);
                    }
                }
            }
        }
        return hsps;
    }

    //1.4
    public static ArrayList<ArrayList<String>> glouton (ArrayList<DBentry> database, String input, ArrayList<ArrayList<Integer>> positionInput, ArrayList<ArrayList<Integer>> positionDatabase, ArrayList<ArrayList<String>> hsp, Integer seuil, String seed){
        ArrayList<ArrayList<String>> hsps = hsp;

        for (int i =0; i < hsps.size(); i++){
            if (!hsp.get(i).isEmpty()) {
                for (int j = 0; j < hsp.get(i).size(); j++) {

                    System.out.println(seed.length());
                    int seqInputFirst = positionInput.get(i).get(j) - 1;
                    int seqInputLast = positionInput.get(i).get(j) + seed.length() + 1;
                    int seqDataFirst = positionDatabase.get(i).get(j)-1;
                    int seqDataLast = positionDatabase.get(i).get(j) + seed.length() + 1;

                    int score = 0;
                    int maxScore = 0;

                    boolean extendLeft = true;
                    boolean extendRight = true;

                    while (extendLeft && extendRight) {
                        if (seqInputFirst == -1 || seqDataFirst == -1 || maxScore - score >= seuil) {
                            extendLeft = false;
                        }

                        if (seqDataLast == database.get(i).getSequence().length() || seqInputLast == input.length() || maxScore - score >= seuil) {
                            extendRight = false;
                        }

                        int scoreLeft = 0;

                        //extend gauche
                        if (extendLeft) {
                            System.out.println(extendLeft);
                            System.out.println(seqInputFirst);
                            System.out.println(seqDataFirst);
                            if (input.charAt(seqInputFirst) == database.get(i).getSequence().charAt(seqDataFirst)){
                                scoreLeft = 5 + score;
                            } else {
                                scoreLeft = -4 + score;
                            }
                        }


                        int scoreRight = 0;

                        //extend droite
                        if (extendRight) {
                            System.out.println(extendRight);
                            System.out.println(seqInputLast);
                            System.out.println(seqDataLast);
                            System.out.println(input);
                            System.out.println(database.get(i).getSequence());
                            if (input.charAt(seqInputLast) == database.get(i).getSequence().charAt(seqDataLast)){
                                scoreRight = 5 + score;
                            } else {
                                scoreRight = -4 + score;
                            }
                        }

                        //on regarde quel est le meilleur extend
                        if (extendRight && scoreLeft >= scoreRight){
                            hsps.get(i).set(j, input.charAt(seqInputFirst) +(hsps.get(i).get(j)));
                            score = scoreLeft;
                            seqInputFirst -= 1;
                            seqDataFirst -= 1;
                        } else if (extendLeft && scoreLeft < scoreRight){
                            hsps.get(i).set(j, (hsps.get(i).get(j))+input.charAt(seqInputLast));
                            score = scoreRight;
                            seqInputLast -= 1;
                            seqDataLast -= 1;
                        }

                        if (score > maxScore){
                            maxScore = score;
                        }
                        System.out.println(hsps);
                    }
                }
            }
        }

        return hsps;
    }
    //1.5
    //1.6
    //1.7
    public static void main(String[] args) throws FileNotFoundException {
        //System.out.println(readFasta("src/tRNAs.fasta"));
        ArrayList<DBentry> db = new ArrayList<DBentry>();
        //ArrayList<DBentry> input = readFasta("src/unknown.fasta");
        DBentry x = new DBentry("BRUH","AAAAB");
        db.add(x);
        DBentry y = new DBentry("BRUH1","BAAAA");
        db.add(y);
        DBentry z = new DBentry("BRUH2","ACAAA");
        db.add(z);
        //readFasta("src/tRNAs.fasta", database);
//        ArrayList<DBentry> unknown = readFasta("src/unknown.fasta");
//        ArrayList<DBentry> database = readFasta("src/tRNAs.fasta");
//        ArrayList<ArrayList<Integer>> position = new ArrayList<ArrayList<Integer>>();
//        ArrayList<ArrayList<Integer>> positionInput = new ArrayList<ArrayList<Integer>>();
//        ArrayList<ArrayList<String>> hsps = hsp(unknown.get(0).getSequence(),database , "11111111111", position, positionInput);
        //System.out.println(kmer(unknown.get(1), 11));

        //Creer ici la database pour avoir acces partout
        //System.out.println("GGGAGAATGACTGAGTGGTTAAAAGTGACAGACTGTAAATCTGTTGAAATTATTTCTACGTAGGTTCGAATCCTGCTTCTCCCA".length());

        //HSP
//        System.out.println(hsp(unknown.get(0).getSequence(),database , "11111111111", position, positionInput));
//
//        //position dans les sequences de la database
//        System.out.println(position);
//
//        //position de l'input
//        System.out.println(positionInput);

        String input = "AAAA";
        ArrayList<ArrayList<Integer>> position = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> positionInput = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<String>> hsps = hsp(input,db,"111", position, positionInput);
        System.out.println(hsps);
        System.out.println(position);
        System.out.println(positionInput);
        System.out.println(glouton(db,input,positionInput, position, hsps,4,"111"));
    }
}
