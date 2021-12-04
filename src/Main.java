import java.io.*;
import java.util.Scanner;
import java.util.ArrayList;

public class Main {
    // 1.1
    //je m'en criss un peu ici so j'ai tout crisser dans un tableau, dans le fond
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
                    boolean bruh = true;
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
                        position.get(y).add(j); // y ieme sequence dans la database , j ieme kmer = index de cette sequence
                        positioninput.get(y).add(i);// y ieme sequence dans la database, i eme kmer = index de notre sequence inconnue
                    }
                }
            }
        }
        return hsps;
    }

    public static void removeDuplicates(ArrayList<ArrayList<String>> list, ArrayList<ArrayList<Integer>> positionInput, ArrayList<ArrayList<Integer>> positionDatabase) {

        // Create a new ArrayList
        ArrayList<ArrayList<String>> newList = new ArrayList<ArrayList<String>>();
        ArrayList<ArrayList<Integer>> positionIn = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> positionDb = new ArrayList<ArrayList<Integer>>();

        // Traverse through the first list
        for (int i = 0;  i < list.size(); i++) {
            ArrayList<String> newist = new ArrayList<String>();
            ArrayList<Integer> posIn = new ArrayList<Integer>();
            ArrayList<Integer> posDb = new ArrayList<Integer>();
            if (!list.get(i).isEmpty()) {
                for (int j = 0; j < list.get(i).size(); j++) {

                    System.out.println(!newist.contains(list.get(i).get(j)) && !posIn.contains(positionInput.get(i).get(j)) && !posDb.contains(positionDatabase.get(i).get(j)));
                    if (!newist.contains(list.get(i).get(j)) || !posIn.contains(positionInput.get(i).get(j)) || !posDb.contains(positionDatabase.get(i).get(j))) {
                        newist.add(list.get(i).get(j));
                        posIn.add(positionInput.get(i).get(j));
                        posDb.add(positionDatabase.get(i).get(j));
                    }
                }
            }
            newList.add(newist);
            positionIn.add(posIn);
            positionDb.add(posDb);
        }

        list.clear();
        positionInput.clear();
        positionDatabase.clear();

        list.addAll(newList);
        positionInput.addAll(positionIn);
        positionDatabase.addAll(positionDb);
    }

    //1.4
    public static ArrayList<ArrayList<String>> glouton (ArrayList<ArrayList<String>> hsp, ArrayList<DBentry> database, String input, ArrayList<ArrayList<Integer>> positionInput, ArrayList<ArrayList<Integer>> positionDatabase, Integer seuil, String seed){

        for (int i = 0; i < hsp.size(); i++){
            if (!hsp.get(i).isEmpty()) {
                for (int j = 0; j < hsp.get(i).size(); j++) {

                    int seqInputFirst = positionInput.get(i).get(j)-1;
                    int seqInputLast = positionInput.get(i).get(j) + seed.length();
                    int seqDataFirst = positionDatabase.get(i).get(j)-1;
                    int seqDataLast = positionDatabase.get(i).get(j) + seed.length();

                    int maxScore = 0;
                    int score = 0;

                    boolean extendLeft = true;
                    boolean extendRight = true;

                    while(extendLeft || extendRight){
                        int scoreLeft = 0;
                        int scoreRight = 0;

                        if (seqInputFirst == -1 || seqDataFirst == -1){
                            extendLeft = false;
                        } else {
                            if (input.charAt(seqInputFirst) == database.get(i).getSequence().charAt(seqDataFirst)){
                                scoreLeft += 5;
                            } else {
                                scoreLeft -= 4;
                            }
                        }

                        if (seqDataLast == database.get(i).getSequence().length()   || seqInputLast == input.length()) {
                            extendRight = false;
                        } else {
                            if (input.charAt(seqInputLast) == database.get(i).getSequence().charAt(seqDataLast)){
                                scoreRight += 5;
                            } else {
                                scoreRight -= 4;
                            }
                        }

                        if (!extendLeft && extendRight) {
                            score += scoreRight;
                            if (maxScore - score < seuil){
                                hsp.get(i).set(j, hsp.get(i).get(j)+input.charAt(seqInputLast));
                                seqInputLast++;
                                seqDataLast++;
                            }
                        }

                        if (extendLeft && !extendRight) {
                            score += scoreLeft;
                            if (maxScore - score < seuil) {
                                hsp.get(i).set(j, input.charAt(seqInputFirst) + hsp.get(i).get(j));
                                positionDatabase.get(i).set(j, positionDatabase.get(i).get(j)-1);
                                positionInput.get(i).set(j, positionInput.get(i).get(j)-1);
                                seqInputFirst--;
                                seqDataFirst--;
                            }
                        }

                        if (extendLeft && extendRight) {
                            score += scoreLeft;
                            if (scoreLeft >= scoreRight){
                                if (maxScore - score < seuil) {
                                    hsp.get(i).set(j, input.charAt(seqInputFirst) + hsp.get(i).get(j));
                                    positionDatabase.get(i).set(j, positionDatabase.get(i).get(j)-1);
                                    positionInput.get(i).set(j, positionInput.get(i).get(j)-1);
                                    seqInputFirst--;
                                    seqDataFirst--;
                                }
                            } else {
                                score += scoreRight;
                                if (maxScore - score < seuil) {
                                    hsp.get(i).set(j, hsp.get(i).get(j) + input.charAt(seqInputLast));
                                    seqInputLast++;
                                    seqDataLast++;
                                }
                            }
                        }

                        if (score > maxScore){
                            maxScore = score;
                        }

                        if (maxScore - score >= seuil){
                            extendLeft = false;
                            extendRight = false;
                        }
                    }
                }
            }
        }

       removeDuplicates(hsp, positionInput, positionDatabase);

        return hsp;
    }
    //1.5
    public static ArrayList<ArrayList<String>> fusionHsp (ArrayList<ArrayList<String>> extendedHsps, ArrayList<ArrayList<Integer>> sequencePositions) {

        //condition de fusion : end index dun hsp > start position dun autre && < end
        for (int i = 0; i < extendedHsps.size(); i++) {
            for (int j = 0; j < extendedHsps.get(i).size(); j++) {

                //sortir si ya aucun hsp a cette ligne ou seulement 1 (donc pas de chevauchement->skip line)
                if(extendedHsps.get(i).size()==0 || extendedHsps.get(i).size()==1){
                    break;
                }

                //end index of current hsp we want to compare
                int hspSize = extendedHsps.get(i).get(j).length();
                int hspEndIndex = sequencePositions.get(i).get(j) + hspSize - 1;

                //chek all other hsps
                for (int l = 0; l < extendedHsps.get(i).size(); l++) {

                    //start and end index of the other hsp
                    int otherHspStart = sequencePositions.get(i).get(l);
                    int otherHspSize = extendedHsps.get(i).get(l).length();
                    int otherHspEnd = sequencePositions.get(i).get(l) + otherHspSize - 1;

                    //condition de fusion
                    if (hspEndIndex > otherHspStart && hspEndIndex < otherHspEnd) {

                        //a partir de end index du current hsp, concat lautre hsp a part sa partie chevauchee
                        int longueurChevauchement = hspEndIndex - otherHspStart;

                        String concatPartOfOtherHsp = (extendedHsps.get(i).get(l).substring(longueurChevauchement));

                        String fusionFinale = extendedHsps.get(i).get(j) + concatPartOfOtherHsp;


                        //add fusion
                        extendedHsps.get(i).set(j, fusionFinale);
                    }
                }
            }
        }
        return extendedHsps;
    }
    //1.6
    //1.7

    public static void plast(String input, String output, String seed, int seuil) throws FileNotFoundException {
        ArrayList<DBentry> unknown = readFasta(input);
        ArrayList<DBentry> database = readFasta(output);
        ArrayList<ArrayList<Integer>> position = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> positionInput = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<String>> hsps = hsp(unknown.get(0).getSequence(), database,seed,position,positionInput);
        System.out.println(database.get(0).getSequence().length());
        ArrayList<ArrayList<String>> hspExtend = glouton(hsps, database, unknown.get(0).getSequence(), positionInput, position, seuil, seed);
        System.out.println(hsps);
        System.out.println(positionInput);
        System.out.println(position);
        //System.out.println(fusionHsp (hsps, position));
    }

    public static void main(String[] args) throws FileNotFoundException {
//        ArrayList<DBentry> db = new ArrayList<DBentry>();
//        DBentry x = new DBentry("BRUH","ABCXXXXABC");
//        db.add(x);
//
//        String input = "PEFABCEFP";
//        ArrayList<ArrayList<Integer>> position = new ArrayList<ArrayList<Integer>>();
//        ArrayList<ArrayList<Integer>> positionInput = new ArrayList<ArrayList<Integer>>();
//        ArrayList<ArrayList<String>> hsps = hsp(input, db,"111",position,positionInput);
//        System.out.println(glouton(hsps, db,input,positionInput, position, 4,"111"));
//        System.out.println(hsps);
//        System.out.println(positionInput);
//        System.out.println(position);

        plast("src/unknown.fasta","src/tRNAs.fasta", "1111", 14);
    }
}
