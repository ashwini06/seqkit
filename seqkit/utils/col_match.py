#!/usr/bin/env python

# merge two files based on the chr-end positions if in two files the 3 cols are chr,start,end
# Usage : colmatch(tss_fl,ndg_fl,cmb_fl,"merge/second_file")

def colmatch(fl1,fl2,fl3,chk):
    pos = {};
    with open(fl1,'r') as in_fl1, open(fl2,'r') as in_fl2, open(fl3, 'w') as op_fl:
        head1 = in_fl1.next(); head2 = in_fl2.next() ;
        head1 = head1.strip().split('\t');
        head2 = head2.strip().split('\t');
        if chk == "merge":
            header = "\t".join(head2+head1)
            op_fl.write(header+"\n")
        elif chk == "second_file" :
            op_fl.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("chr","start","end","length","abs_summit","pileup","-log10(pvalue)","fold_enrichment","-log10(qvalue)","name"))
        for l in in_fl1:
            ln = l.strip().split('\t')
            # store chr and end pos
            id1 = "-".join([ln[0],ln[2]])
            pos[id1] = ln
        for m in in_fl2:
            mn = m.strip()
            if re.search(r'^chr.*',mn):
                mn = mn.split('\t')
                id2 = "-".join([mn[0],mn[2]])
                if id2 in pos :
                    if chk == "merge" :
                        op = "\t".join(mn+pos[id2])
                        op_fl.write(op+"\n")

