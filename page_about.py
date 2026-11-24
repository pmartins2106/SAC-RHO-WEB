
import streamlit as st

def about():
    
    # st.markdown("<h2 style='text-align: center;'>About</h2>", 
    #             unsafe_allow_html=True)
    
    st.subheader("📖 Guidelines2")
    st.markdown(
        "**SAC-RHO** extracts the **Seeding Activity Coefficient (SAC)** and seed concentration parameter "
        "**RHO** ($\\rho$) from ThT fluorescence aggregation curves automatically. Main steps:"
    )
    # st.write('**SAC-RHO** is a free-to-use webtool for the analysis of protein aggregation data such as Thioflavin-T (ThT) fluorescence traces. Primarily oriented for seed amplification assays like the RT-QuIC assay, the principles of SAC-RHO (available [here]( http://www.jbc.org/lookup/doi/10.1074/jbc.M112.375345), [here](https://onlinelibrary.wiley.com/doi/10.1002/advs.202301501) and [here](https://www.molbiolcell.org/doi/10.1091/mbc.E23-07-0289)) also apply to any phase separation process dominated by secondary nucleation and/or growth. Raw data describing mass-based progress curves are used as input to provide seeding activity coefficients (SACs) and RHO estimates as output. RHO values linearly increase with seed concentration.')
    
    # Lower next markdowns
    st.write("")
    
    with st.expander("📊 1. Data Preparation & Upload"):
        st.markdown(
            """
            * **Format:** Save data as an OpenDocument Spreadsheet (`.ods`).
            * **Structure:** Place **Time** in the first column, followed by your data columns.
            * **Flexibility:** Varying column lengths and manual cell deletion (to remove outliers) are fully supported.
            """
        )

    with st.expander("⚙️ 2. Configuration Settings (Defaults Optimized)"):
        st.markdown(
            """
            * **Nucleation Rate ($k_b$):** Default is `0` (seeded growth). Increase if spontaneous nucleation occurs.
            * **Preestablished SAC:** Leave blank to calculate from the upload, or input a fixed baseline value to compare datasets.
            * **Reference Endpoint ($F_\\infty^0$):** Automatically uses your highest plateau. Input a fixed value to standardize across different datasets.
            * **Minimum Datapoints ($N_{min}$):** Default is `15`. Ensures enough points are analyzed in the early growth phase (initial 20%).
            """
        )

    with st.expander("📈 3. Output & Results"):
        st.markdown(
            """
            * **Visual Fits:** Plots raw data alongside fitted trendlines to visually confirm accuracy.
            * **Downloadable Table:** Contains calculated **SAC** (dataset baseline), **RHO** (relative seed concentrations), and quality metrics with standard errors to spot outliers.
            """
        )
    
    def make_line():
        """ Line divider between images. """
            
        line = st.markdown('<hr style="border:1px solid gray"> </hr>',
                unsafe_allow_html=True)

        return line    
   
    # make_line()
    
    # Lower next markdowns
    # st.write("")

    
    # st.markdown("<h2 style='text-align: center;'>Disclaimer</h2>", 
    #             unsafe_allow_html=True)
    st.subheader("⚠️ Disclaimer")
    st.markdown('''
    This webtool is free to use for all non-commercial purposes. This webtool is provided by the copyright holders and contributors "as is" and any express or implied warranties, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose are disclaimed. In no event shall the copyright owner or contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) however caused and on any theory of liability, whether in contract, strict liability, or tort (including negligence or otherwise) arising in any way out of the use of this webtool, even if advised of the possibility of such damage.
    ''')
    
    make_line()
    
    # Lower next markdowns
    st.write("")
    
    # st.markdown("<h2 style='text-align: center;'>Help and Support</h2>", 
    #             unsafe_allow_html=True)
    st.subheader("💬 Help and Support")
    st.markdown('''
    If you need technical assistance or a stand-alone version of ACCU-RATES please contact
    pmartins@i3s.up.pt or frocha@fe.up.pt. 
    ''')
    
    make_line()
    
    # Lower next markdowns
    st.write("")
    
    # st.markdown("<h2 style='text-align: center;'>Funding</h2>", 
    #             unsafe_allow_html=True)
    
    # st.image("images//logos.png", use_container_width=True)
    st.subheader("🤝 Funding")
    st.markdown('''
                This work is funded by the Portuguese Foundation for Science and Technology (FCT) in the framework of projects FCT 2023.11892.PEX and COMPETE2030-FEDER-00690800.
    ''')


       

    
