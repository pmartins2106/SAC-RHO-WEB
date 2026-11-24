"""
SAC-RHO: Webtool for automatic analysis of seed amplification assays.
"""

from typing import Callable, Dict
import streamlit as st

from page_about import about
from page_SACRHO import analyse


def main() -> None:
    """Configure the Streamlit sidebar and handle page navigation."""
    st.set_page_config(page_title="SAC-RHO", page_icon="🧬")

    st.sidebar.success(
        "🚀 Webtool for automatic analysis of seed amplification assays"
    )
    st.sidebar.title("Main options")

    pages: Dict[str, Callable[[], None]] = {
        "ANALYSE": analyse,
        "GUIDELINES": about,
    }

    selected_page = st.sidebar.radio("Select:", list(pages.keys()))
    pages[selected_page]()


if __name__ == "__main__":
    main()
