import sys
import argparse
import ast
    
def get_dt_token(dt: float) -> str:
    """
    Converts a float dt into a string token suitable for filenames.
    E.g., 1e-6 -> "1e6", 2.5e-7 -> "25e8"
    """
    mantissa, exponent = f"{dt:e}".split("e")
    mantissa = mantissa.rstrip("0").rstrip(".") or "0"
    exponent = exponent.lstrip("+-")
    exponent = exponent.lstrip("0") or "0"
    dt_token = f"{mantissa}e{exponent}"
    # print(f"DT token: {dt_token}")
    return dt_token

def get_viscosity_token(viscosity: float) -> str:
    """
    Converts a float viscosity into a string token suitable for filenames.
    E.g., 0.003 -> "003", 0.025 -> "025"
    """
    token = int(round(viscosity * 1e9))
    token_str = f"{token:09d}"
    # remove trailing zeros
    token_str = token_str.rstrip("0")
    # print(f"Viscosity token: {token_str}")
    return token_str

def main():
    available_functions = {
        "main": main,
        "get_viscosity_token": get_viscosity_token,
        "get_dt_token": get_dt_token,
    }

    parser = argparse.ArgumentParser(description="Execute functions from main.py")
    parser.add_argument("func", nargs="?", default="main", choices=available_functions.keys())
    parser.add_argument("func_args", nargs="*", help="Positional or key=value arguments")
    parsed = parser.parse_args()

    def _coerce(token: str):
        try:
            return ast.literal_eval(token)
        except (ValueError, SyntaxError):
            return token

    positional, keyword = [], {}
    for token in parsed.func_args:
        if "=" in token:
            key, value = token.split("=", 1)
            keyword[key] = _coerce(value)
        else:
            positional.append(_coerce(token))

    try:
        available_functions[parsed.func](*positional, **keyword)
    except KeyboardInterrupt:
        print("\nProgram interrupted by user. Exiting gracefully.")
        sys.exit(1)

if __name__ == "__main__":
    main()