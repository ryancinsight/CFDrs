import os
import toml

def check_workspace():
    if not os.path.exists("Cargo.toml"):
        print("Root Cargo.toml not found")
        return

    try:
        data = toml.load("Cargo.toml")
    except Exception as e:
        print(f"Failed to parse Cargo.toml: {e}")
        return

    members = data.get("workspace", {}).get("members", [])
    print(f"Checking {len(members)} members...")

    for member in members:
        member_path = os.path.join(os.getcwd(), member)
        cargo_path = os.path.join(member_path, "Cargo.toml")
        if not os.path.exists(cargo_path):
            print(f"MISSING: {member} -> {cargo_path}")
        else:
            print(f"OK: {member}")

if __name__ == "__main__":
    check_workspace()
