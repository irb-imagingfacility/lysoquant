import re
def tokenize(s):
    return re.findall(r'"[^"]*"|[{}:]|[^\s{}:]+', s)
def parse_block(tok, i):
    d = {}
    while i < len(tok) and tok[i] != '}':
        k = tok[i]; i += 1
        if tok[i] == ':':
            v = tok[i+1]; i += 2
            v = v.strip('"')
            d.setdefault(k, []).append(v)
        else:
            assert tok[i] == '{', (k, tok[i])
            sub, i = parse_block(tok, i+1)
            i += 1
            d.setdefault(k, []).append(sub)
    return d, i
def load(path):
    tok = tokenize(open(path).read())
    d, _ = parse_block(tok, 0)
    return d
