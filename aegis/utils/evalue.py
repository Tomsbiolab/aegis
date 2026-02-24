def parse_evalue(e):
    e = e.strip()
    try:
        if e.startswith('>'):
            return float(round_evalue(e[1:]))  # handle cases like '>10'
        elif e.lower() in ('na', 'nan', ''):
            return float('nan')  # handle missing values
        else:
            return float(round_evalue(e))
            
    except ValueError:
        print(f"Error: Could not parse E-value: {e}")
        return float('nan')  # or raise error, depending on use case

def round_evalue(e):
    e = float(e)
    e = f"{e:.2e}"

    return e