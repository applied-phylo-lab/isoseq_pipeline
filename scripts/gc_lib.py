"""
Shared helpers for the gene conversion analysis stage.

The V gene FASTA headers written by make_vgene_fasta.py encode the locus
coordinates we need for the topology model:

    {prefix}_{locus}.{pos}.{contig}.{gene_type}.{productive}.{strand}

e.g. VGP_redwinged_blackbird_IGL.8749176.CM036732.1.V.True.-
"""
import re
from collections import namedtuple

# ─── FASTA ────────────────────────────────────────────────────────────────────

def read_fasta(path):
    """Return dict: id -> sequence (uppercase). Ids are truncated at whitespace."""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif name is not None:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n{seq}\n")


# ─── V gene names ─────────────────────────────────────────────────────────────

GeneInfo = namedtuple("GeneInfo", "name locus pos contig productive strand")

# ...IGL.8749176.CM036732.1.V.True.-   (contig itself contains dots)
_GENE_RE = re.compile(
    r"^(?P<prefix>.+)_(?P<locus>IG[HKL]|TR[ABGD])\."
    r"(?P<pos>\d+)\.(?P<contig>.+)\.(?P<type>[VDJ])\."
    r"(?P<productive>True|False)\.(?P<strand>[+-])(?:_\d+)?$"
)


def parse_gene_name(name):
    """Parse a pipeline V gene FASTA header into a GeneInfo, or None."""
    m = _GENE_RE.match(name)
    if not m:
        return None
    return GeneInfo(
        name=name,
        locus=m.group("locus"),
        pos=int(m.group("pos")),
        contig=m.group("contig"),
        productive=(m.group("productive") == "True"),
        strand=m.group("strand"),
    )


def parse_gene_names(names):
    """Parse many names; raises if any fail so topology is never silently wrong."""
    out, bad = {}, []
    for n in names:
        g = parse_gene_name(n)
        if g is None:
            bad.append(n)
        else:
            out[n] = g
    if bad:
        raise SystemExit(
            "Could not parse locus coordinates from these V gene names:\n  "
            + "\n  ".join(bad[:10])
            + ("\n  ..." if len(bad) > 10 else "")
        )
    return out


# ─── PAF + cs tag ─────────────────────────────────────────────────────────────

PafRec = namedtuple(
    "PafRec",
    "query qlen qstart qend strand target tlen tstart tend nmatch alnlen mapq cs",
)


def parse_paf(path):
    """Yield PafRec for each alignment line, with the cs tag if present."""
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            cs = None
            for tag in f[12:]:
                if tag.startswith("cs:Z:"):
                    cs = tag[5:]
                    break
            yield PafRec(
                query=f[0], qlen=int(f[1]), qstart=int(f[2]), qend=int(f[3]),
                strand=f[4], target=f[5], tlen=int(f[6]),
                tstart=int(f[7]), tend=int(f[8]),
                nmatch=int(f[9]), alnlen=int(f[10]), mapq=int(f[11]), cs=cs,
            )


_CS_RE = re.compile(r"(:\d+|\*[acgtn]{2}|\+[acgtn]+|-[acgtn]+|=[ACGTN]+)")


def cs_to_pairs(rec):
    """
    Walk a minimap2 cs tag and return a list of (target_pos, query_base) for
    every aligned target position, target_pos being 0-based on the V gene.

    Substitutions give the query base; matches give the target base; positions
    deleted from the target are reported as '-'. Insertions are skipped (they
    do not occupy a target coordinate).

    Returns None when the record has no cs tag.
    """
    if rec.cs is None:
        return None
    tpos = rec.tstart
    pairs = []
    for tok in _CS_RE.findall(rec.cs):
        op = tok[0]
        if op == ":":                      # run of identity, bases not spelled out
            n = int(tok[1:])
            pairs.extend((tpos + i, None) for i in range(n))
            tpos += n
        elif op == "=":                    # run of identity, bases spelled out
            for i, b in enumerate(tok[1:]):
                pairs.append((tpos + i, b.upper()))
            tpos += len(tok) - 1
        elif op == "*":                    # substitution: ref, then query
            pairs.append((tpos, tok[2].upper()))
            tpos += 1
        elif op == "-":                    # deletion from the target
            for i in range(len(tok) - 1):
                pairs.append((tpos + i, "-"))
            tpos += len(tok) - 1
        elif op == "+":                    # insertion, no target coordinate
            pass
    return pairs


def projected_query(rec, target_seq):
    """
    Project an alignment onto V gene coordinates.

    Returns a list the length of the V gene, holding the transcript's base at
    each V gene position, or None where the alignment does not cover it.
    Identity runs in the cs tag are filled from target_seq.
    """
    pairs = cs_to_pairs(rec)
    if pairs is None:
        return None
    out = [None] * rec.tlen
    for tpos, qbase in pairs:
        if 0 <= tpos < rec.tlen:
            out[tpos] = target_seq[tpos] if qbase is None else qbase
    return out
