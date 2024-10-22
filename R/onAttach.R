.onAttach <- function(...) {
    msg <- paste(
        "  Dusa, Adrian (2019) QCA with R. A Comprehensive Resource.",
        "\n",
        "  Springer International Publishing.",
        "\n\n",
        "To run the graphical user interface, use: runGUI()",
        "\n",
        sep = ""
    )

    packageStartupMessage("\nTo cite package QCA in publications, please use:\n", msg, "\n")
}

