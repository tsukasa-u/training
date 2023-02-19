#![allow(dead_code, unused_imports, unused_variables)]


fn type_of<T>(_: T) -> String{
    let a = std::any::type_name::<T>();
    return a.to_string();
}


extern crate proc_macro;
use quote::quote;
use syn::{ext::IdentExt, parse_macro_input};

#[proc_macro_attribute]
pub fn test1(attr: proc_macro::TokenStream, item: proc_macro::TokenStream) -> proc_macro::TokenStream {
    println!("{}", attr);
    println!("{:?}", attr);
    println!("{}", item);
    println!("{:?}", item);
    return proc_macro::TokenStream::default();
}

#[derive(Debug)]
struct TestMacroInput {
    from: syn::LitInt,
    to: syn::LitInt,
    inclusive: bool,
    ident: syn::Ident,
    tt: proc_macro2::TokenStream,
}

impl syn::parse::Parse for TestMacroInput {
    fn parse(input: syn::parse::ParseStream) -> syn::Result<Self> {
        println!("{:?}", input);
        let ident = syn::Ident::parse(input)?;
        let _in = <syn::Token![in]>::parse(input)?;
        // let ch = syn::ExprBinary::parse(input)?;
        let from = syn::LitInt::parse(input)?;
        let inclusive = input.peek(syn::Token![..=]);
        if inclusive {
            <syn::Token![..=]>::parse(input)?;
        } else {
            <syn::Token![..]>::parse(input)?;
        }
        let to = syn::LitInt::parse(input)?;
        let content;
        let _braces = syn::braced!(content in input);
        let tt = proc_macro2::TokenStream::parse(&content)?;

        Ok(TestMacroInput {
            from,
            to,
            inclusive,
            ident,
            tt
        })
    }
}

#[proc_macro]
pub fn test2(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    // Parse the input tokens into a syntax tree
    let input = parse_macro_input!(input as TestMacroInput);
    println!("{:?}", input);

    // Build the output, possibly using quasi-quotation
    let expanded = quote! {
        // ...
    };

    // Hand the output tokens back to the compiler
    return proc_macro::TokenStream::from(expanded);
}

